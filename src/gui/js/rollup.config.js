import resolve from '@rollup/plugin-node-resolve';
import commonjs from '@rollup/plugin-commonjs';
import terser from '@rollup/plugin-terser';

export default {
  input: 'src/index.js',
  output: {
    file: 'dist/bundle.js',
    format: 'esm'
  },
  plugins: [
    resolve(),
    commonjs(),
    terser({
      compress: {
        drop_console: true,
        drop_debugger: true
      },
      mangle: true,
      format: {
        beautify: true
      }
    })
  ],
  external: ['three'],
  treeshake: true
};