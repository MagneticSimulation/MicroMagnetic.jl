import resolve from '@rollup/plugin-node-resolve';
import commonjs from '@rollup/plugin-commonjs';
import terser from '@rollup/plugin-terser';

export default {
  input: 'src/index.js',
  output: [
    {
      file: 'dist/bundle.js',
      format: 'umd',
      name: 'MicromagneticGUI',
      sourcemap: false
    },
    {
      file: 'dist/bundle.esm.js',
      format: 'esm',
      sourcemap: false
    }
  ],
  treeshake: true,
  plugins: [
    resolve(),
    commonjs(),
    terser({
      compress: {
        drop_console: true,
        drop_debugger: true
      },
      mangle: true,
      output: {
        beautify: true,
        comments: false
      }
    })
  ]
};