// PrusaSlicer is released under the terms of the AGPLv3 or higher.

#include "TransferExporter.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <nlohmann/json.hpp>
#include <tiffio.h>

#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/SLA/RasterBase.hpp"
#include "libslic3r/libslic3r.h"

namespace Slic3r {
namespace transfer {
namespace {

constexpr size_t max_dimension_px = 32768;
// The rasterizer uses one byte per pixel, before packing the TIFF bits.
constexpr size_t max_raster_pixels = 100 * 1024 * 1024;

size_t raster_dimension(double length_mm, double pixel_size_mm)
{
    const double pixels = length_mm / pixel_size_mm;
    if (!std::isfinite(pixels) || pixels <= 0 || pixels > max_dimension_px)
        throw std::runtime_error("Transfer raster dimensions must be between 1 and 32768 pixels.");

    // Avoid an extra pixel caused only by floating-point rounding.
    const double nearest = std::round(pixels);
    return static_cast<size_t>(nearest >= 1 && std::abs(pixels - nearest) < 1e-9 ?
                               nearest : std::ceil(pixels));
}

void write_bilevel_tiff(const void *data, size_t width, size_t height,
                       size_t components, const TransferExportConfig &config,
                       const std::string &filename,
                       const std::string &description)
{
    if (data == nullptr || components != 1 || width != config.width_px ||
        height != config.height_px)
        throw std::runtime_error("Unexpected transfer raster layout.");

    // B explicitly selects MSB-first bits in the scanline buffers.
    std::unique_ptr<TIFF, decltype(&TIFFClose)> tif(TIFFOpen(filename.c_str(), "wB"), &TIFFClose);
    if (!tif)
        throw std::runtime_error("Cannot open transfer TIFF: " + filename);

    if (TIFFSetField(tif.get(), TIFFTAG_IMAGEWIDTH, static_cast<uint32_t>(width)) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_IMAGELENGTH, static_cast<uint32_t>(height)) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_BITSPERSAMPLE, 1) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_SAMPLESPERPIXEL, 1) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_COMPRESSION, COMPRESSION_NONE) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_FILLORDER, FILLORDER_MSB2LSB) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_ROWSPERSTRIP, static_cast<uint32_t>(64)) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_RESOLUTIONUNIT, RESUNIT_CENTIMETER) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_XRESOLUTION, 10.0 / config.pixel_size_mm) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_YRESOLUTION, 10.0 / config.pixel_size_mm) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_IMAGEDESCRIPTION, description.c_str()) != 1 ||
        TIFFSetField(tif.get(), TIFFTAG_SOFTWARE, "PrusaSlicer Transfer Export") != 1)
        throw std::runtime_error("Cannot set transfer TIFF fields: " + filename);

    const auto *pixels = static_cast<const uint8_t *>(data);
    std::vector<uint8_t> scanline((width + 7) / 8, 0);
    for (size_t y = 0; y < height; ++y) {
        std::fill(scanline.begin(), scanline.end(), uint8_t(0));
        for (size_t x = 0; x < width; ++x) {
            if (pixels[y * width + x] >= 128)
                scanline[x / 8] |= static_cast<uint8_t>(0x80u >> (x % 8));
        }
        if (TIFFWriteScanline(tif.get(), scanline.data(), static_cast<uint32_t>(y), 0) != 1)
            throw std::runtime_error("Cannot write transfer TIFF scanline: " + filename);
    }

    if (TIFFWriteDirectory(tif.get()) != 1 || TIFFFlush(tif.get()) != 1)
        throw std::runtime_error("Cannot finalize transfer TIFF: " + filename);
}

} // namespace

TransferExportConfig load_transfer_export_config(const Points &bed_shape,
                                                const std::string &config_path)
{
    std::ifstream input(config_path);
    if (!input)
        throw std::runtime_error("Cannot read transfer export configuration: " + config_path);

    nlohmann::json json;
    try {
        json = nlohmann::json::parse(input);
    } catch (const nlohmann::json::exception &error) {
        throw std::runtime_error("Invalid transfer export configuration " + config_path +
                                 ": " + error.what());
    }
    if (!json.is_object() || json.size() != 1 || !json.contains("pixel_size_mm") ||
        !json.at("pixel_size_mm").is_number())
        throw std::runtime_error("Transfer configuration must contain only numeric pixel_size_mm.");

    const double pixel_size_mm = json.at("pixel_size_mm").get<double>();
    if (!std::isfinite(pixel_size_mm) || pixel_size_mm <= 0)
        throw std::runtime_error("Transfer pixel_size_mm must be finite and greater than zero.");

    const BoundingBox plate_bounds(bed_shape);
    if (bed_shape.size() != 4 || !plate_bounds.defined)
        throw std::runtime_error("Transfer export requires an axis-aligned rectangular bed.");

    for (size_t i = 0; i < bed_shape.size(); ++i) {
        const Point &point = bed_shape[i];
        const Point &next = bed_shape[(i + 1) % bed_shape.size()];
        const Point &opposite = bed_shape[(i + 2) % bed_shape.size()];
        const bool corner = (point.x() == plate_bounds.min.x() ||
                             point.x() == plate_bounds.max.x()) &&
                            (point.y() == plate_bounds.min.y() ||
                             point.y() == plate_bounds.max.y());
        const bool axis_aligned_edge = (point.x() == next.x()) != (point.y() == next.y());
        if (!corner || !axis_aligned_edge || point == opposite)
            throw std::runtime_error("Transfer export requires an axis-aligned rectangular bed.");
    }

    const Point plate_size = plate_bounds.max - plate_bounds.min;
    const size_t width_px = raster_dimension(unscale<double>(plate_size.x()), pixel_size_mm);
    const size_t height_px = raster_dimension(unscale<double>(plate_size.y()), pixel_size_mm);
    if (width_px > max_raster_pixels / height_px)
        throw std::runtime_error("Transfer raster exceeds the 100 MiB pixel buffer limit.");

    return {plate_bounds, pixel_size_mm, width_px, height_px};
}

std::string create_transfer_export_directory()
{
    const boost::filesystem::path base_dir = boost::filesystem::absolute(
        boost::filesystem::path("output") / "transfer_print");
    boost::filesystem::create_directories(base_dir);

    const std::tm local_time = boost::posix_time::to_tm(
        boost::posix_time::second_clock::local_time());
    std::ostringstream job_name;
    job_name << "job_" << std::put_time(&local_time, "%Y-%m-%d_%H%M%S");

    for (size_t suffix = 0; ; ++suffix) {
        std::ostringstream name;
        name << job_name.str();
        if (suffix > 0)
            name << '_' << std::setw(2) << std::setfill('0') << suffix;

        const boost::filesystem::path job_dir = base_dir / name.str();
        // Directory creation reserves the name, including for concurrent exports.
        if (boost::filesystem::create_directory(job_dir))
            return job_dir.string();
    }
}

std::string export_transfer_layer(const Layer &layer, const Point &instance_shift,
                                 const TransferExportConfig &config,
                                 const std::string &output_directory)
{
    // Work on copies; slice coordinates are centered around the PrintObject.
    ExPolygons polygons = layer.lslices;
    translate(polygons, instance_shift);

    for (const ExPolygon &polygon : polygons) {
        if (polygon.contour.points.size() < 3)
            throw std::runtime_error("Invalid transfer contour.");
        for (const Polygon &hole : polygon.holes)
            if (hole.points.size() < 3)
                throw std::runtime_error("Invalid transfer hole.");
        if (!config.plate_bounds.contains(get_extents(polygon)))
            throw std::runtime_error(
                "Layer is outside the transfer plate defined by the printer profile.");
    }

    // The image starts at the plate's minimum X/Y, which need not be zero.
    translate(polygons, Point(-config.plate_bounds.min.x(), -config.plate_bounds.min.y()));

    // The default transform maps plate +Y upwards to image rows downwards.
    // Gamma zero thresholds the raster coverage, producing binary pixel values.
    auto raster = sla::create_raster_grayscale_aa(
        sla::Resolution(config.width_px, config.height_px),
        sla::PixelDim(config.pixel_size_mm, config.pixel_size_mm), 0.0);
    for (const ExPolygon &polygon : polygons)
        raster->draw(polygon);

    std::ostringstream name;
    name << "layer_" << std::setw(4) << std::setfill('0') << layer.id() << ".tif";
    const std::string filename =
        (boost::filesystem::path(output_directory) / name.str()).string();

    const Point plate_size = config.plate_bounds.max - config.plate_bounds.min;
    std::ostringstream description;
    description << std::setprecision(12) << "Transfer export; plate_mm="
                << unscale<double>(plate_size.x()) << 'x'
                << unscale<double>(plate_size.y())
                << "; image_px=" << config.width_px << 'x' << config.height_px
                << "; pixel_size_mm=" << config.pixel_size_mm
                << "; origin_mm=" << unscale<double>(config.plate_bounds.min.x()) << ','
                << unscale<double>(config.plate_bounds.min.y())
                << "; plate_axes=+X right,+Y up; white=apply; black=background"
                << "; layer_index=" << layer.id()
                << "; slice_z_mm=" << layer.slice_z
                << "; print_z_mm=" << layer.print_z
                << "; layer_height_mm=" << layer.height;
    const std::string metadata = description.str();

    // The encoder callback exposes the raw raster buffer without an image file.
    // TIFF is written directly to disk, so no encoded memory buffer is returned.
    raster->encode([&](const void *data, size_t width, size_t height, size_t components) {
        write_bilevel_tiff(data, width, height, components, config, filename, metadata);
        return sla::EncodedRaster{};
    });

    return filename;
}

} // namespace transfer
} // namespace Slic3r
