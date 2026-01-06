/**
 * Colormap utilities for scientific visualization
 * Provides common colormaps used in research plots for surface coloring
 */

/**
 * Clamp a value between min and max
 */
function clamp(value, min = 0, max = 1) {
    return Math.max(min, Math.min(max, value));
}

/**
 * Convert HSV to RGB color
 * @param {number} h - Hue (0-1)
 * @param {number} s - Saturation (0-1)
 * @param {number} v - Value/Brightness (0-1)
 * @returns {Object} RGB color {r, g, b} in range 0-1
 */
function hsvToRgb(h, s, v) {
    let r, g, b;
    const i = Math.floor(h * 6);
    const f = h * 6 - i;
    const p = v * (1 - s);
    const q = v * (1 - f * s);
    const t = v * (1 - (1 - f) * s);

    switch (i % 6) {
        case 0: r = v; g = t; b = p; break;
        case 1: r = q; g = v; b = p; break;
        case 2: r = p; g = v; b = t; break;
        case 3: r = p; g = q; b = v; break;
        case 4: r = t; g = p; b = v; break;
        case 5: r = v; g = p; b = q; break;
    }

    return { r, g, b };
}

/**
 * Coolwarm colormap (blue-white-red)
 * Blue for negative, white for zero, red for positive values
 * Good for showing divergence around zero
 */
export const coolwarm = {
    name: 'Coolwarm',
    description: 'Blue-white-red divergent colormap',
    getColor(t) {
        t = clamp(t);
        // Blue (0) -> White (0.5) -> Red (1)
        if (t < 0.5) {
            // Blue to white
            const s = t * 2;
            return {
                r: s,
                g: s,
                b: 1
            };
        } else {
            // White to red
            const s = (t - 0.5) * 2;
            return {
                r: 1,
                g: 1 - s * 0.5,
                b: 1 - s * 0.5
            };
        }
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Blue-White-Red colormap (same as coolwarm)
 */
export const bwr = coolwarm;

/**
 * Red-White-Blue divergent colormap (reversed coolwarm)
 */
export const rwb = {
    name: 'Red-White-Blue',
    description: 'Red-white-blue divergent colormap',
    getColor(t) {
        t = clamp(t);
        if (t < 0.5) {
            const s = t * 2;
            return {
                r: 1,
                g: 1 - s * 0.5,
                b: 1 - s * 0.5
            };
        } else {
            const s = (t - 0.5) * 2;
            return {
                r: 1 - s,
                g: 1 - s,
                b: 1
            };
        }
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * HSV/Rainbow colormap
 * Traditional rainbow colormap with full hue spectrum
 */
export const hsv = {
    name: 'HSV (Rainbow)',
    description: 'Full hue spectrum from HSV color space',
    getColor(t) {
        t = clamp(t);
        // Use full hue range from 0.66 (blue) to 0.0 (red) for better visibility
        // But standard HSV uses 0-1 mapping to hue circle
        return hsvToRgb(1 - t, 1, 1);
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Jet colormap
 * Classic colormap used in many visualization tools
 * Dark blue -> blue -> cyan -> green -> yellow -> red
 */
export const jet = {
    name: 'Jet',
    description: 'Classic jet colormap (dark blue to red)',
    getColor(t) {
        t = clamp(t);
        const n = 7;
        const i = Math.floor(t * (n - 1));
        const f = t * (n - 1) - i;

        const colors = [
            { r: 0, g: 0, b: 0.5 },      // Dark blue
            { r: 0, g: 0, b: 1 },        // Blue
            { r: 0, g: 1, b: 1 },        // Cyan
            { r: 0, g: 1, b: 0 },        // Green
            { r: 1, g: 1, b: 0 },        // Yellow
            { r: 1, g: 0, b: 0 },        // Red
            { r: 0.5, g: 0, b: 0 }       // Dark red
        ];

        const c1 = colors[i];
        const c2 = colors[Math.min(i + 1, n - 1)];

        return {
            r: c1.r + (c2.r - c1.r) * f,
            g: c1.g + (c2.g - c1.g) * f,
            b: c1.b + (c2.b - c1.b) * f
        };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Viridis colormap
 * Perceptually uniform colormap, colorblind-friendly
 * Dark blue -> purple -> green -> yellow
 */
export const viridis = {
    name: 'Viridis',
    description: 'Perceptually uniform, colorblind-friendly colormap',
    getColor(t) {
        t = clamp(t);
        // Viridis RGB values at specific points
        const colors = [
            { r: 0.267004, g: 0.004874, b: 0.329415 },
            { r: 0.282327, g: 0.140926, b: 0.457517 },
            { r: 0.253935, g: 0.265254, b: 0.529983 },
            { r: 0.206756, g: 0.371758, b: 0.553117 },
            { r: 0.163625, g: 0.466133, b: 0.550347 },
            { r: 0.127568, g: 0.566949, b: 0.550556 },
            { r: 0.134692, g: 0.658636, b: 0.517649 },
            { r: 0.266941, g: 0.748751, b: 0.440573 },
            { r: 0.477504, g: 0.821444, b: 0.318195 },
            { r: 0.741388, g: 0.873449, b: 0.149561 },
            { r: 0.993248, g: 0.906157, b: 0.143936 }
        ];

        const n = colors.length - 1;
        const i = Math.floor(t * n);
        const f = t * n - i;

        if (i >= n) return colors[n];

        const c1 = colors[i];
        const c2 = colors[i + 1];

        return {
            r: c1.r + (c2.r - c1.r) * f,
            g: c1.g + (c2.g - c1.g) * f,
            b: c1.b + (c2.b - c1.b) * f
        };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Plasma colormap
 * Similar to viridis, high contrast
 * Dark blue -> purple -> red -> yellow
 */
export const plasma = {
    name: 'Plasma',
    description: 'High contrast, perceptually uniform colormap',
    getColor(t) {
        t = clamp(t);
        const colors = [
            { r: 0.050383, g: 0.029803, b: 0.527976 },
            { r: 0.089219, g: 0.064884, b: 0.586004 },
            { r: 0.129783, g: 0.110705, b: 0.638209 },
            { r: 0.178909, g: 0.169335, b: 0.685151 },
            { r: 0.240346, g: 0.231171, b: 0.720434 },
            { r: 0.307398, g: 0.295733, b: 0.746186 },
            { r: 0.380095, g: 0.360693, b: 0.758802 },
            { r: 0.460316, g: 0.428399, b: 0.754698 },
            { r: 0.546448, g: 0.498062, b: 0.734991 },
            { r: 0.633964, g: 0.566623, b: 0.701627 },
            { r: 0.716613, g: 0.633164, b: 0.657958 },
            { r: 0.792418, g: 0.695619, b: 0.603644 },
            { r: 0.859324, g: 0.752527, b: 0.538602 },
            { r: 0.915080, g: 0.801332, b: 0.464022 },
            { r: 0.957610, g: 0.840605, b: 0.381094 },
            { r: 0.982541, g: 0.870427, b: 0.291005 },
            { r: 0.988883, g: 0.891705, b: 0.196198 }
        ];

        const n = colors.length - 1;
        const i = Math.floor(t * n);
        const f = t * n - i;

        if (i >= n) return colors[n];

        const c1 = colors[i];
        const c2 = colors[i + 1];

        return {
            r: c1.r + (c2.r - c1.r) * f,
            g: c1.g + (c2.g - c1.g) * f,
            b: c1.b + (c2.b - c1.b) * f
        };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Inferno colormap
 * Dark through orange to yellow
 * Good for dark backgrounds
 */
export const inferno = {
    name: 'Inferno',
    description: 'Dark to bright colormap, good for dark backgrounds',
    getColor(t) {
        t = clamp(t);
        const colors = [
            { r: 0.001462, g: 0.000466, b: 0.013866 },
            { r: 0.014462, g: 0.005032, b: 0.067406 },
            { r: 0.038756, g: 0.016926, b: 0.133888 },
            { r: 0.074425, g: 0.036608, b: 0.209132 },
            { r: 0.119702, g: 0.064244, b: 0.284536 },
            { r: 0.172857, g: 0.098315, b: 0.355844 },
            { r: 0.232163, g: 0.138116, b: 0.419552 },
            { r: 0.295798, g: 0.183776, b: 0.473088 },
            { r: 0.362326, g: 0.234511, b: 0.514037 },
            { r: 0.430458, g: 0.288735, b: 0.540068 },
            { r: 0.498863, g: 0.344234, b: 0.549862 },
            { r: 0.565981, g: 0.398706, b: 0.543088 },
            { r: 0.628728, g: 0.450874, b: 0.522201 },
            { r: 0.684539, g: 0.499004, b: 0.489878 },
            { r: 0.732327, g: 0.541894, b: 0.454189 },
            { r: 0.770455, g: 0.578722, b: 0.420383 },
            { r: 0.798983, g: 0.608351, b: 0.390183 },
            { r: 0.818558, g: 0.632021, b: 0.365947 },
            { r: 0.836879, g: 0.654539, b: 0.347988 },
            { r: 0.859023, g: 0.679974, b: 0.334529 },
            { r: 0.885038, g: 0.709681, b: 0.326148 },
            { r: 0.910641, g: 0.740208, b: 0.325102 },
            { r: 0.931836, g: 0.769893, b: 0.332055 },
            { r: 0.949105, g: 0.798717, b: 0.347300 },
            { r: 0.962642, g: 0.825866, b: 0.370243 },
            { r: 0.973962, g: 0.851320, b: 0.400028 },
            { r: 0.983246, g: 0.874836, b: 0.435655 },
            { r: 0.990306, g: 0.896605, b: 0.477103 },
            { r: 0.995605, g: 0.916229, b: 0.521014 },
            { r: 0.998770, g: 0.933594, b: 0.566156 },
            { r: 1.000000, g: 0.948000, b: 0.608000 }
        ];

        const n = colors.length - 1;
        const i = Math.floor(t * n);
        const f = t * n - i;

        if (i >= n) return colors[n];

        const c1 = colors[i];
        const c2 = colors[i + 1];

        return {
            r: c1.r + (c2.r - c1.r) * f,
            g: c1.g + (c2.g - c1.g) * f,
            b: c1.b + (c2.b - c1.b) * f
        };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Gray colormap
 * Grayscale from black to white
 */
export const gray = {
    name: 'Gray',
    description: 'Grayscale from black to white',
    getColor(t) {
        t = clamp(t);
        return { r: t, g: t, b: t };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};


/**
 * Hot colormap
 * Black -> red -> yellow -> white
 */
export const hot = {
    name: 'Hot',
    description: 'Black-red-yellow-white colormap',
    getColor(t) {
        t = clamp(t);
        if (t < 1/3) {
            return { r: t * 3, g: 0, b: 0 };
        } else if (t < 2/3) {
            return { r: 1, g: (t - 1/3) * 3, b: 0 };
        } else {
            return { r: 1, g: 1, b: (t - 2/3) * 3 };
        }
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Cool colormap
 * Cyan-magenta gradient
 */
export const cool = {
    name: 'Cool',
    description: 'Cyan-magenta gradient',
    getColor(t) {
        t = clamp(t);
        return { r: t, g: 1 - t, b: 1 };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Spring colormap
 * Magenta-yellow gradient
 */
export const spring = {
    name: 'Spring',
    description: 'Magenta-yellow gradient',
    getColor(t) {
        t = clamp(t);
        return { r: 1, g: t, b: 1 - t };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Autumn colormap
 * Red-yellow gradient
 */
export const autumn = {
    name: 'Autumn',
    description: 'Red-yellow gradient',
    getColor(t) {
        t = clamp(t);
        return { r: 1, g: t, b: 0 };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Winter colormap
 * Blue-green gradient
 */
export const winter = {
    name: 'Winter',
    description: 'Blue-green gradient',
    getColor(t) {
        t = clamp(t);
        return { r: 0, g: t, b: 1 - t };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Copper colormap
 * Black-copper gradient
 */
export const copper = {
    name: 'Copper',
    description: 'Black-copper gradient',
    getColor(t) {
        t = clamp(t);
        const f = Math.sqrt(t);
        return { r: f * 0.9 + 0.1 * t, g: f * 0.5, b: f * 0.3 };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Bone colormap
 * Black-white with blue tint (like X-ray)
 */
export const bone = {
    name: 'Bone',
    description: 'Black-white with blue tint',
    getColor(t) {
        t = clamp(t);
        return {
            r: t * 0.9,
            g: t * 0.9,
            b: t * 1.0
        };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Spectral colormap
 * Blue-green-yellow-red (perceptually uniform)
 */
export const spectral = {
    name: 'Spectral',
    description: 'Blue-green-yellow-red perceptually uniform colormap',
    getColor(t) {
        t = clamp(t);
        const colors = [
            { r: 0.158, g: 0.317, b: 0.584 },
            { r: 0.368, g: 0.501, b: 0.729 },
            { r: 0.554, g: 0.681, b: 0.812 },
            { r: 0.742, g: 0.836, b: 0.655 },
            { r: 0.993, g: 0.906, b: 0.144 }
        ];

        const n = colors.length - 1;
        const i = Math.floor(t * n);
        const f = t * n - i;

        if (i >= n) return colors[n];

        const c1 = colors[i];
        const c2 = colors[i + 1];

        return {
            r: c1.r + (c2.r - c1.r) * f,
            g: c1.g + (c2.g - c1.g) * f,
            b: c1.b + (c2.b - c1.b) * f
        };
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Seismic colormap
 * Blue-white-red (similar to coolwarm but more saturated)
 */
export const seismic = {
    name: 'Seismic',
    description: 'Blue-white-red for seismic data',
    getColor(t) {
        t = clamp(t);
        if (t < 0.5) {
            const s = t * 2;
            return {
                r: s * 0.5,
                g: s * 0.5,
                b: s * 0.8 + (1 - s) * 0.2
            };
        } else {
            const s = (t - 0.5) * 2;
            return {
                r: s * 0.8 + (1 - s) * 0.5,
                g: s * 0.2 + (1 - s) * 0.5,
                b: (1 - s) * 0.8
            };
        }
    },
    getColorRGB(t) {
        const c = this.getColor(t);
        return {
            r: Math.round(c.r * 255),
            g: Math.round(c.g * 255),
            b: Math.round(c.b * 255)
        };
    }
};

/**
 * Get array of colors from a colormap
 * @param {Object} colormap - Colormap object
 * @param {number} numColors - Number of colors to generate
 * @returns {Array} Array of {r, g, b} color objects
 */
export function getColorArray(colormap, numColors = 256) {
    const colors = [];
    for (let i = 0; i < numColors; i++) {
        colors.push(colormap.getColor(i / (numColors - 1)));
    }
    return colors;
}

/**
 * Get RGB array for use in WebGL
 * @param {Object} colormap - Colormap object
 * @param {number} numColors - Number of colors
 * @returns {Float32Array} Flat array of RGB values
 */
export function getColorArrayRGB(colormap, numColors = 256) {
    const colors = new Float32Array(numColors * 3);
    for (let i = 0; i < numColors; i++) {
        const c = colormap.getColor(i / (numColors - 1));
        colors[i * 3] = c.r;
        colors[i * 3 + 1] = c.g;
        colors[i * 3 + 2] = c.b;
    }
    return colors;
}

/**
 * Get a color from a value using a specified colormap
 * @param {number} value - Value to map (0-1)
 * @param {string|Object} colormapName - Name of colormap or colormap object
 * @returns {Object} {r, g, b} color in range 0-1
 */
export function getColor(value, colormapName = 'coolwarm') {
    const colormap = typeof colormapName === 'string' 
        ? colormaps[colormapName.toLowerCase()] 
        : colormapName;
    
    if (!colormap) {
        console.warn(`Unknown colormap: ${colormapName}, using coolwarm`);
        return coolwarm.getColor(clamp(value));
    }
    
    return colormap.getColor(clamp(value));
}

/**
 * Normalize a value to 0-1 range
 * @param {number} value - Value to normalize
 * @param {number} min - Minimum value
 * @param {number} max - Maximum value
 * @returns {number} Normalized value (0-1)
 */
export function normalizeValue(value, min, max) {
    if (max === min) return 0.5;
    return clamp((value - min) / (max - min));
}

/**
 * Map values to colors using a colormap
 * @param {Array} values - Array of values to map
 * @param {string|Object} colormapName - Name of colormap
 * @param {number} min - Minimum value (optional, auto-calculated if not provided)
 * @param {number} max - Maximum value (optional, auto-calculated if not provided)
 * @returns {Array} Array of {r, g, b} colors
 */
export function mapValuesToColors(values, colormapName = 'coolwarm', min, max) {
    const actualMin = min !== undefined ? min : Math.min(...values);
    const actualMax = max !== undefined ? max : Math.max(...values);
    
    return values.map(v => getColor(normalizeValue(v, actualMin, actualMax), colormapName));
}

/**
 * Export all available colormaps
 */
export const colormaps = {
    coolwarm,
    bwr,
    rwb,
    hsv,
    jet,
    viridis,
    plasma,
    inferno,
    gray,
    hot,
    cool,
    spring,
    autumn,
    winter,
    copper,
    bone,
    spectral,
    seismic
};

/**
 * Get list of available colormap names
 */
export function getAvailableColormaps() {
    return Object.keys(colormaps);
}

/**
 * Get a specific colormap by name
 */
export function getColormap(name) {
    return colormaps[name.toLowerCase()] || coolwarm;
}

export default {
    coolwarm,
    bwr,
    rwb,
    hsv,
    jet,
    viridis,
    plasma,
    inferno,
    gray,
    hot,
    cool,
    spring,
    autumn,
    winter,
    copper,
    bone,
    spectral,
    seismic,
    colormaps,
    getAvailableColormaps,
    getColormap,
    getColor,
    getColorArray,
    getColorArrayRGB,
    normalizeValue,
    mapValuesToColors,
    hsvToRgb
};
