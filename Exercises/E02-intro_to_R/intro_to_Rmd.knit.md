---
title: "E02 - Intro to Rmd (Analysis of Gene Expression @ VSCHT)"
author:
  - Jiri Novotny jiri.novotny@img.cas.cz
  - Studuj bioinformatiku! http://studuj.bioinformatiku.cz
institute: "Laboratory of Genomics and Bioinformatics @ Institute of Molecular Genetics of the ASCR"
output:
  html_document:
    toc: true
    toc_depth: 4
    toc_float: true
date: "2021-03-09"
---

***



# RMarkdown (Rmd)

We will be using [RMarkdown](https://rmarkdown.rstudio.com/lesson-1.html) for exercises.
It allows to combine [Markdown](https://pandoc.org/MANUAL.html#pandocs-markdown) and (not only) R code
to generate documents in HTML, PDF, Word (docx) and many more formats.
This is [how it works](https://rmarkdown.rstudio.com/lesson-2.html):

![Code chunks in RMarkdown document are evaluated by [knitr](https://yihui.name/knitr/) and converted to [Markdown](https://pandoc.org/MANUAL.html#pandocs-markdown),
which is then compiled by [pandoc](https://pandoc.org) to selected output format (in our case HTML).
All of this is handled automatically by `rmarkdown::render()` function (see the end of this document).](/data/persistent/jirinovo/bio-class-deb10/AGE2021/Exercises/E02-intro_to_R/images/rmarkdown_flow.png)

**Useful links**

- [How RMarkdown works](https://uclouvain-cbio.github.io/WSBIM1207/sec-rr.html)
- [R Markdown: The Definitive Guide](https://bookdown.org/yihui/rmarkdown/)
- [RMarkdown cheatsheet](https://rmarkdown.rstudio.com/lesson-15.html)
  (PDF [here](https://raw.githubusercontent.com/rstudio/cheatsheets/master/rmarkdown-2.0.pdf) and
  [here](https://rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf))
- [RMarkdown cookbook](https://bookdown.org/yihui/rmarkdown-cookbook/)
- [Markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
- [pandoc's Markdown flavour](https://pandoc.org/MANUAL.html#pandocs-markdown)
- [knitr options](https://yihui.name/knitr/options/)
- [RMarkdown output formats](https://rmarkdown.rstudio.com/lesson-9.html)

***

## How to use RMarkdown

**Pandoc header**

On the beginning of the source file you can see the header section separated by `---`.
Those are options for `pandoc`, mostly specifying how output should look like.
They are written in [YAML format](https://lzone.de/cheat-sheet/YAML).

You can see we are using three common parameters: `title`, `author` or `institute`.
The first two will be shown in the header of resulting document.

As the output type we are using `html_document` with several parameters:

- `toc: true`: make a Table of Contents.
- `toc_depth: 4`: how many section levels to show in the TOC.
- `toc_float: true`: make TOC floating on the left side of page. If `false`, TOC will be placed statically to the beginning of document.

For possible output types see [this](https://bookdown.org/yihui/rmarkdown/documents.html) chapter from R Markdown: The Definitive Guide.

It is also possible to use inline R code in parameters. As an example, we used this to output current date:

<pre class="r">date: "&grave;r Sys.Date()&grave;"</pre>

**Markdown syntax**

Markdown is an easy-to-write/easy-to-read markup language (see the source of this document).
It is widely used by, but not only, developers for documentation (you have certainly seen it on GitHub).
Markdown syntax is super easy and so, please, refer to [Markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
(however, `pandoc` uses an [extended version](https://pandoc.org/MANUAL.html#pandocs-markdown) of Markdown).

**Using HTML tags**

Besides Markdown syntax, you can use anything you would normally use in HTML.
For example, `<div style="color: red;">Some red text.</div>`:

<div style="color: red;">Some red text.</div>

Just be careful, because RMarkdown is not a classic HTML document and you don't have a full control over conversion to output HTML.

**Code chunks**

Code (R, Python, etc.) is placed in [chunks](https://bookdown.org/yihui/rmarkdown/r-code.html).
A chunk begins with triple backticks followed by `{language_name}` and ends with another triple backticks.
Both start and end triple backticks must be on separate lines.

Example of chunk with R code:

<pre class="r">
&grave;&grave;&grave;{r}
df <- data.frame()
&grave;&grave;&grave;
</pre>

Output from a chunk is generated below the chunk.

Each chunk can have options set. For example,

<pre class="r">
&grave;&grave;&grave;{r, warning = FALSE}
df <- data.frame()
&grave;&grave;&grave;
</pre>

will suppress all warnings produced by this chunk. See [knitr options](https://yihui.name/knitr/options/).

You can set these options globally via `knitr` package:


```r
library(knitr)
opts_chunk$set(warning = FALSE)
```

And override them in individual chunks.

JavaScript and CSS are also [supported](https://bookdown.org/yihui/rmarkdown/language-engines.html#javascript-and-css).
Code in chunks of type `js` and `css` will be put to `<script>` and `<style>` tags inside `<head>`.

**Inline code**

You can also place R code [inline](https://rmarkdown.rstudio.com/lesson-4.html) using <code  class="r">``` `r code` ```</code>.
Inline code in final document will be replaced by its output.

**Notebook mode**

In Rstudio, you use RMarkdown documents as so-called [notebooks](https://bookdown.org/yihui/rmarkdown/notebook.html).
That is, each code chunk can be evaluated in current R session. Let's see some examples (run them in your RStudio):


```r
head(mtcars)
```


```r
plot(mpg ~ cyl, data = mtcars)
```

<img src="/data/persistent/jirinovo/bio-class-deb10/AGE2021/Exercises/E02-intro_to_R/intro_to_Rmd_files/figure-html/unnamed-chunk-5-1.png" width="672" />

**Rendering a document**

On the end of this document you will find a chunk, which generates a HTML document from this source.

## Useful Rmd-related keyboard shortcuts in RStudio

You can view them in `Tools -> Keyboard Shortcuts Help`.

- `Ctrl + Enter`: run current line or selected lines.
- `Ctrl + Alt + R`: run all chunks.
- `Ctrl + Alt + P`: run all chunks above.
- `Ctrl + Shift + Enter`: run current chunk.
- `Ctrl + Alt + N`: run next chunk.

## Some RMarkdown tips

### Pretty tables

We can use [DT](https://rstudio.github.io/DT/) package, an R interface to JavaScript library [datatables](https://datatables.net/).


```r
DT::datatable(
  mtcars[, 1:5],
  filter = "top",
  class = "display",
  rownames = TRUE
)
```

```{=html}
<div id="htmlwidget-a0f33e4ab27445afecd0" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-a0f33e4ab27445afecd0">{"x":{"filter":"top","filterHTML":"<tr>\n  <td><\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"10.4\" data-max=\"33.9\" data-scale=\"1\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"4\" data-max=\"8\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"71.1\" data-max=\"472\" data-scale=\"1\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"52\" data-max=\"335\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"2.76\" data-max=\"4.93\" data-scale=\"2\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","data":[["Mazda RX4","Mazda RX4 Wag","Datsun 710","Hornet 4 Drive","Hornet Sportabout","Valiant","Duster 360","Merc 240D","Merc 230","Merc 280","Merc 280C","Merc 450SE","Merc 450SL","Merc 450SLC","Cadillac Fleetwood","Lincoln Continental","Chrysler Imperial","Fiat 128","Honda Civic","Toyota Corolla","Toyota Corona","Dodge Challenger","AMC Javelin","Camaro Z28","Pontiac Firebird","Fiat X1-9","Porsche 914-2","Lotus Europa","Ford Pantera L","Ferrari Dino","Maserati Bora","Volvo 142E"],[21,21,22.8,21.4,18.7,18.1,14.3,24.4,22.8,19.2,17.8,16.4,17.3,15.2,10.4,10.4,14.7,32.4,30.4,33.9,21.5,15.5,15.2,13.3,19.2,27.3,26,30.4,15.8,19.7,15,21.4],[6,6,4,6,8,6,8,4,4,6,6,8,8,8,8,8,8,4,4,4,4,8,8,8,8,4,4,4,8,6,8,4],[160,160,108,258,360,225,360,146.7,140.8,167.6,167.6,275.8,275.8,275.8,472,460,440,78.7,75.7,71.1,120.1,318,304,350,400,79,120.3,95.1,351,145,301,121],[110,110,93,110,175,105,245,62,95,123,123,180,180,180,205,215,230,66,52,65,97,150,150,245,175,66,91,113,264,175,335,109],[3.9,3.9,3.85,3.08,3.15,2.76,3.21,3.69,3.92,3.92,3.92,3.07,3.07,3.07,2.93,3,3.23,4.08,4.93,4.22,3.7,2.76,3.15,3.73,3.08,4.08,4.43,3.77,4.22,3.62,3.54,4.11]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>mpg<\/th>\n      <th>cyl<\/th>\n      <th>disp<\/th>\n      <th>hp<\/th>\n      <th>drat<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":[],"jsHooks":[]}</script>
```

Or [kable](https://bookdown.org/yihui/rmarkdown-cookbook/kable.html) and
[kableExtra](https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html) for simplier tables.


```r
kableExtra::kable_styling(
  kableExtra::kbl(head(mtcars)),
  bootstrap_options = c("striped", "hover", "condensed")
)
```

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> mpg </th>
   <th style="text-align:right;"> cyl </th>
   <th style="text-align:right;"> disp </th>
   <th style="text-align:right;"> hp </th>
   <th style="text-align:right;"> drat </th>
   <th style="text-align:right;"> wt </th>
   <th style="text-align:right;"> qsec </th>
   <th style="text-align:right;"> vs </th>
   <th style="text-align:right;"> am </th>
   <th style="text-align:right;"> gear </th>
   <th style="text-align:right;"> carb </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Mazda RX4 </td>
   <td style="text-align:right;"> 21.0 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 160 </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 3.90 </td>
   <td style="text-align:right;"> 2.620 </td>
   <td style="text-align:right;"> 16.46 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mazda RX4 Wag </td>
   <td style="text-align:right;"> 21.0 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 160 </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 3.90 </td>
   <td style="text-align:right;"> 2.875 </td>
   <td style="text-align:right;"> 17.02 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Datsun 710 </td>
   <td style="text-align:right;"> 22.8 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 108 </td>
   <td style="text-align:right;"> 93 </td>
   <td style="text-align:right;"> 3.85 </td>
   <td style="text-align:right;"> 2.320 </td>
   <td style="text-align:right;"> 18.61 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hornet 4 Drive </td>
   <td style="text-align:right;"> 21.4 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 258 </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 3.08 </td>
   <td style="text-align:right;"> 3.215 </td>
   <td style="text-align:right;"> 19.44 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hornet Sportabout </td>
   <td style="text-align:right;"> 18.7 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 360 </td>
   <td style="text-align:right;"> 175 </td>
   <td style="text-align:right;"> 3.15 </td>
   <td style="text-align:right;"> 3.440 </td>
   <td style="text-align:right;"> 17.02 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Valiant </td>
   <td style="text-align:right;"> 18.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 225 </td>
   <td style="text-align:right;"> 105 </td>
   <td style="text-align:right;"> 2.76 </td>
   <td style="text-align:right;"> 3.460 </td>
   <td style="text-align:right;"> 20.22 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>

***

### Using tab design

As you can see, this document has autogenerated TOC menu on the left.
However, you can also use a tab-style menu.
See the `tab_design.Rmd` ([HTML](tab_design.html)) example.

Personally, I have found this design very useful for reporting.

***

### Generating Markdown and HTML from R

The cool feature of code chunks is also the ability of **raw** output, by setting `results = "asis"`.
For example, this chunk produces a Markdown header and some text below:

<pre class="r">
&grave;&grave;&grave;{r, results = "asis"}
cat("\n\n#### I am an example header\n")
cat("\nHello world!\n\n")
&grave;&grave;&grave;
</pre>


#### I am an example header

Hello world!

I have found this to be very useful for dynamic reports, e.g. header names are composed from sample groups.

***

# Rendering RMarkdown

To render a Rmd file, use the `rmarkdown::render()` function.
Below you can find a code chunk, which loads the `knitr` and `here` packages, sets global chunk options and renders a HTML file.


```r
library(knitr)
library(here)

# You can set global chunk options. Options set in individual chunks will override this.
opts_chunk$set(warning = FALSE, message = FALSE, eval = TRUE)
rmarkdown::render(
  here("E02-intro_to_R/intro_to_Rmd.Rmd"),
  output_file = here("E02-intro_to_R/intro_to_Rmd.html"),
  clean = FALSE
)
```

Because `clean = FALSE` in `rmarkdown::render()`, intermediate Markdown (`.md`) files are not removed and you have an opportunity to inspect them.
These files are sources for `pandoc`, which does the actual rendering to HTML format (as specified in the header of this Rmd file, but
`pandoc` can render to many other formats than HTML, e.g. PDF or LaTex).

## Working directory in a Rmd file

Unfortunately, there are discrepancies in the current working directory when you evaluate chunks in RStudio vs. render a Rmd file,
because `knitr` will automatically set the working directory to that of Rmd file location.

When you evaluate this chunk in RStudio, you get the working directory of your current project.


```r
current_wd <- getwd()
current_wd
```

```
## [1] "/data/persistent/jirinovo/bio-class-deb10/AGE2021/Exercises/E02-intro_to_R"
```

But as you can see in the output `intro_to_Rmd.html`, the working directory is `Exercises/E02-intro_to_R` there.

So how to use consistent project-relative paths in a Rmd file?
Use the [here](https://github.com/jennybc/here_here) package, which is described in `intro_to_R.Rmd` in Reproducible R section.