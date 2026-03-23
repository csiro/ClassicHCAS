// Verify that float is 32-bit ensuring accuracy and vectorised perfromance
// Don't want to rely on C++23 just for float32_t
static_assert(sizeof(float) == 4, "float must be 32-bit");
static_assert(std::numeric_limits<float>::is_iec559, "float must be IEEE 754");
// Create float32_t for clarity
using float32_t = float;
