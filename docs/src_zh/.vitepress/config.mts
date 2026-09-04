import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import { withMermaid } from "vitepress-plugin-mermaid";

import path from 'path'
import fs from 'fs'

const mathjax = mathjaxPlugin()

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

const nav = [
  ...navTemp.nav,
  {
    text: 'English',
    link: 'https://magneticsimulation.github.io/MicroMagnetic.jl/dev/'
  }
]

// https://vitepress.dev/reference/site-config
export default withMermaid (defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  // Docstrings (English, shared with the English site) contain section @refs like
  // `[Micromagnetic model](@ref)` that only resolve against English page headings.
  // Ignore just those unresolved `@ref` hrefs; genuine dead links still fail the build.
  ignoreDeadLinks: [/\/@ref/],
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
