import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import { withMermaid } from "vitepress-plugin-mermaid";

import path from 'path'
import fs from 'node:fs'

const mathjax = mathjaxPlugin()

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
    link: `${getBaseRepository(baseTemp.base)}dev/`
  }
]

// docs/TODO.md #3 (Option B): before the Vitepress bundling step, create placeholder
// files for any referenced `public/` media that was not generated at build time
// (e.g. an example block failed or ran in draft mode), so that a single page cannot
// break the whole build with "Could not resolve". The make script exports
// DOCUMENTER_MD_ROOT (the markdown output directory, e.g. `build_zh/.documenter`).
const mediaPlaceholders: Record<string, Buffer> = {
  png: Buffer.from('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg==', 'base64'),
  jpg: Buffer.from('/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAgGBgcGBQgHBwcJCQgKDBQNDAsLDBkSEw8UHRofHh0aHBwgJC4nICIsIxwcKDcpLDAxNDQ0Hyc5PTgyPDUzNDP/wAALCAABAAEBAREA/8QAFAABAAAAAAAAAAAAAAAAAAAACf/EABQQAQAAAAAAAAAAAAAAAAAAAAD/2gAIAQEAAD8AKp//2Q==', 'base64'),
  jpeg: Buffer.from('/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAgGBgcGBQgHBwcJCQgKDBQNDAsLDBkSEw8UHRofHh0aHBwgJC4nICIsIxwcKDcpLDAxNDQ0Hyc5PTgyPDUzNDP/wAALCAABAAEBAREA/8QAFAABAAAAAAAAAAAAAAAAAAAACf/EABQQAQAAAAAAAAAAAAAAAAAAAAD/2gAIAQEAAD8AKp//2Q==', 'base64'),
  gif: Buffer.from('R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7', 'base64'),
  svg: Buffer.from('<svg xmlns="http://www.w3.org/2000/svg" width="8" height="8"><rect width="8" height="8" fill="#ddd"/></svg>'),
  pdf: Buffer.from('%PDF-1.4\n')
}

function fillMissingMedia(): void {
  const mdRoot = process.env.DOCUMENTER_MD_ROOT
  if (!mdRoot || !fs.existsSync(mdRoot)) return
  const exts = 'mp4|webm|png|jpe?g|gif|svg|pdf'
  const publicRegex = new RegExp(`public\\/[\\w.\\-/]+\\.(?:${exts})`, 'g')
  const imgRegex = new RegExp(`!\\[[^\\]]*\\]\\(([^)\\s]+\\.(?:${exts}))\\)`, 'g')
  const imgHtmlRegex = new RegExp(`<img[^>]+src=["']([^"'\s]+\\.(?:${exts}))["']`, 'g')
  const targets = new Set<string>()
  const walk = (dir: string) => {
    for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
      if (e.name === '.vitepress' || e.name === 'public' || e.name === 'node_modules') continue
      const p = path.join(dir, e.name)
      if (e.isDirectory()) walk(p)
      else if (e.name.endsWith('.md')) {
        const text = fs.readFileSync(p, 'utf8')
        for (const m of text.matchAll(publicRegex)) targets.add(path.join(mdRoot, 'public', m[0].slice('public/'.length)))
        for (const m of text.matchAll(imgRegex)) {
          if (m[1].startsWith('http') || m[1].startsWith('/')) continue
          targets.add(path.resolve(dir, m[1]))
        }
        let rewritten = text
        for (const m of text.matchAll(imgHtmlRegex)) {
          const src = m[1]
          if (src.startsWith('http') || src.startsWith('/') || src.startsWith('./') || src.startsWith('../')) continue
          targets.add(path.resolve(dir, src))
          // Draft-mode Documenter emits raw `<img src="x.png">`; a bare specifier is
          // resolved as a module id by Vite and fails even when the file exists.
          rewritten = rewritten.split(`src="${src}"`).join(`src="./${src}"`)
          rewritten = rewritten.split(`src='${src}'`).join(`src='./${src}'`)
        }
        if (rewritten !== text) fs.writeFileSync(p, rewritten)
      }
    }
  }
  walk(mdRoot)
  const missing: string[] = []
  for (const target of targets) {
    if (fs.existsSync(target)) continue
    fs.mkdirSync(path.dirname(target), { recursive: true })
    const ext = target.split('.').pop()!.toLowerCase()
    fs.writeFileSync(target, mediaPlaceholders[ext] ?? Buffer.alloc(0))
    missing.push(path.relative(mdRoot, target))
  }
  if (missing.length > 0) {
    console.warn(`[fillMissingMedia] created ${missing.length} placeholder file(s) for missing media: ${missing.join(', ')}`)
  }
}
fillMissingMedia()

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
