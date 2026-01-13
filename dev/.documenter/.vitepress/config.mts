import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// const nav = [
//   ...navTemp.nav,
//   {
//     component: 'VersionPicker'
//   }
// ]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/GeneFinder.jl/dev/',// TODO: replace this in makedocs!
  title: 'GeneFinder.jl',
  description: "A VitePress Site",
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }]],
  ignoreDeadLinks: true,

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    logo: { src: '/logo.svg', width: 24, height: 24},
    // logo: { src: '/logo.png', width: 24, height: 24 },
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav: [
{ text: 'Home', link: '/index' },
{ text: 'Get Started', link: '/getstarted' },
{ text: 'Usage', collapsed: false, items: [
{ text: 'The ORF type', link: '/orftype' },
{ text: 'Scoring ORFs', link: '/features' },
{ text: 'A Simple Coding Rule', link: '/simplecodingrule' },
{ text: 'Ribosome Binding Sites', link: '/rbs' },
{ text: 'Writing ORFs In Files', link: '/iodocs' }]
 },
{ text: 'API', link: '/api' }
]
,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Get Started', link: '/getstarted' },
{ text: 'Usage', collapsed: false, items: [
{ text: 'The ORF type', link: '/orftype' },
{ text: 'Scoring ORFs', link: '/features' },
{ text: 'A Simple Coding Rule', link: '/simplecodingrule' },
{ text: 'Ribosome Binding Sites', link: '/rbs' },
{ text: 'Writing ORFs In Files', link: '/iodocs' }]
 },
{ text: 'API', link: '/api' }
]
,
    editLink: { pattern: "https://https://github.com/camilogarciabotero/GeneFinder.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/camilogarciabotero/GeneFinder.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    }
  }
})