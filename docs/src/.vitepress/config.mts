import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import { withMermaid } from "vitepress-plugin-mermaid";

import path from 'path'
import fs from 'fs'

const mathjax = mathjaxPlugin()

// Draft builds (DOCS_DRAFT=true) skip executing example blocks, so media that
// examples generate at build time (e.g. eigen_fig2.png, public/*.mp4) is
// missing and Vite's rollup fails on every unresolved import. Scan the
// generated markdown for local media references and materialise placeholder
// files; also normalise the bare <img src="x.png"> paths Documenter emits in
// draft mode to ./x.png. In full builds the files exist and nothing is written.
const PNG_1X1 = Buffer.from(
  'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==',
  'base64')

function fillMissingMedia(root: string) {
  const walk = (dir: string, acc: string[] = []): string[] => {
    for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
      const p = path.join(dir, e.name)
      if (e.isDirectory()) {
        if (e.name !== 'node_modules' && e.name !== '.vitepress') walk(p, acc)
      } else if (e.name.endsWith('.md')) {
        acc.push(p)
      }
    }
    return acc
  }
  const placeholder = (p: string, kind: 'img' | 'video') => {
    if (fs.existsSync(p)) return
    fs.mkdirSync(path.dirname(p), { recursive: true })
    fs.writeFileSync(p, kind === 'img' ? PNG_1X1 : Buffer.alloc(0))
    console.warn(`[fillMissingMedia] placeholder: ${path.relative(root, p)}`)
  }
  const isLocal = (ref: string) => !/^[a-z]+:\/\//.test(ref) && !ref.startsWith('/') && !ref.startsWith('#')
  for (const f of walk(root)) {
    let txt = fs.readFileSync(f, 'utf8')
    const orig = txt
    txt = txt.replace(/(<(?:img|video|source)\s[^>]*?src=")([^"\/:][^"]*)(")/g,
                      (m, a, src, z) => `${a}./${src}${z}`)
    for (const m of txt.matchAll(/!\[[^\]]*\]\(([^)#?]+)\)/g)) {
      const ref = m[1].trim()
      if (!isLocal(ref)) continue
      const target = path.resolve(path.dirname(f), decodeURIComponent(ref))
      if (/\.(png|jpe?g|gif|svg)$/i.test(target)) placeholder(target, 'img')
      else if (/\.(mp4|webm)$/i.test(target)) placeholder(target, 'video')
    }
    for (const m of txt.matchAll(/<img\s[^>]*?src="([^"]+)"/g)) {
      const ref = m[1].trim()
      if (!isLocal(ref)) continue
      const target = path.resolve(path.dirname(f), decodeURIComponent(ref))
      if (/\.(png|jpe?g|gif|svg)$/i.test(target)) placeholder(target, 'img')
    }
    for (const m of txt.matchAll(/<(?:video|source)\s[^>]*?src="([^"]+)"/g)) {
      const ref = m[1].trim()
      if (!isLocal(ref)) continue
      const target = path.resolve(path.dirname(f), decodeURIComponent(ref))
      if (/\.(mp4|webm)$/i.test(target)) placeholder(target, 'video')
    }
    if (txt !== orig) fs.writeFileSync(f, txt)
  }
}

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

// The generated navbar lists every top-level section; fold the example
// sections into one dropdown to keep the bar narrow.
const navGenerated = [...navTemp.nav]
const exampleSections = ['Atomistic', 'Micromagnetics (FD)', 'Micromagnetics (FE)', 'Miscellaneous']
const nav = [
  ...(() => {
    const items = navGenerated.filter((e) => exampleSections.includes(e.text)).flatMap((e) => e.items ?? [])
    const folded = [...navGenerated.filter((e) => !exampleSections.includes(e.text))]
    let i = folded.findIndex((e) => e.text === 'User Guide')
    if (i < 0) i = folded.findIndex((e) => e.text === 'API') - 1
    folded.splice(i + 1, 0, { text: 'Examples', items })
    return folded
  })(),
  {
    text: '中文',
    link: 'https://magneticsimulation.github.io/MicroMagnetic.jl/zh/'
  },
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default withMermaid (defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}],
  ],

  vite: {
    plugins: [
      mathjax.vitePlugin,
      {
        name: 'fill-missing-media',
        buildStart() {
          fillMissingMedia(path.resolve(__dirname, '..'))
        },
      },
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ],
    },
    ssr: {
      noExternal: [
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ],
    },
  },
  markdown: {
    codeTransformers: [juliaReplTransformer()],
    config(md) {
      md.use(tabsMarkdownPlugin);
      md.use(footnote);
      mathjax.markdownConfig(md);
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },

  themeConfig: {
    outline: 'deep',
    logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    },
  }
}))
