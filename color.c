#include "color.h"

color_t colorFromUint32(uint32_t c) {
    return (color_t) {
        (uint8_t) ((c & 0x00FF0000) >> 16),
        (uint8_t) ((c & 0x0000FF00) >> 8),
        (uint8_t) (c & 0x000000FF)
    };
}

color_t sumColors3(color_t c0, color_t c1, color_t c2) {
    return (color_t) {
        (uint8_t) clamp(c0.r + c1.r + c2.r, 255.0),
        (uint8_t) clamp(c0.g + c1.g + c2.g, 255.0),
        (uint8_t) clamp(c0.b + c1.b + c2.b, 255.0)
    };
}

color_t sumColors(color_t c0, color_t c1) {
    return (color_t) {
        (uint8_t) clamp(c0.r + c1.r, 255.0),
        (uint8_t) clamp(c0.g + c1.g, 255.0),
        (uint8_t) clamp(c0.b + c1.b, 255.0)
    };
}

color_t colorFromFloats(float r, float g, float b) {
    return (color_t) {
        (uint8_t) clamp(r * 255.0, 255.0),
        (uint8_t) clamp(g * 255.0, 255.0),
        (uint8_t) clamp(b * 255.0, 255.0)
    };
}