// PrusaSlicer is released under the terms of the AGPLv3 or higher.

#ifndef slic3r_TransferExporter_hpp_
#define slic3r_TransferExporter_hpp_

#include <cstddef>
#include <string>

#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"

namespace Slic3r {

class Layer;

namespace transfer {

struct TransferExportConfig {
    BoundingBox plate_bounds;
    double pixel_size_mm;
    size_t width_px;
    size_t height_px;
};

// Read pixel_size_mm from JSON and derive the raster from a rectangular bed.
// Dimensions round up to cover the plate without changing the pixel size.
TransferExportConfig load_transfer_export_config(const Points &bed_shape,
                                                const std::string &config_path);

// Create a job directory named job_YYYY-MM-DD_HHMMSS using local time below
// output/transfer_print in the process working directory. Existing names get
// a numeric suffix (_01, _02, ...). The returned path is absolute.
std::string create_transfer_export_directory();

// Export one layer using the plate bounds and raster dimensions in config.
// instance_shift is in PrusaSlicer's scaled coordinates. output_directory must
// already exist; all layers of a job share this directory.
std::string export_transfer_layer(const Layer &layer, const Point &instance_shift,
                                 const TransferExportConfig &config,
                                 const std::string &output_directory);

} // namespace transfer
} // namespace Slic3r

#endif
