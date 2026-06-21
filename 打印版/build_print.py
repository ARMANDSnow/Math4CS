# -*- coding: utf-8 -*-
"""Markdown -> 紧密打印 HTML (KaTeX 渲染数学)。自己从磁盘读源文件。"""
import re, html, sys

SRC = {
    "小抄": "/Users/dingyuxuan/Desktop/学校课程/Math4CS/全量小抄.md",
    "作业讲解": "/Users/dingyuxuan/Desktop/学校课程/Math4CS/作业全解-HW1-6.md",
}
OUT = {
    "小抄": "/Users/dingyuxuan/Desktop/小抄-打印版.html",
    "作业讲解": "/Users/dingyuxuan/Desktop/作业讲解-打印版.html",
}
FONT = {"小抄": "7.7px", "作业讲解": "5.3px"}
COLS = {"小抄": "3", "作业讲解": "3"}
LH   = {"小抄": "1.18", "作业讲解": "1.12"}
MARGIN = {"小抄": "5mm", "作业讲解": "4mm"}
# 按文档注入的额外样式:作业讲解去掉每题 h2 的黑块,改为黑字下划线
_NOBAR = ("h2{{background:none!important;color:#000!important;padding:1.5px 0 0.5px!important;"
          "border:0!important;border-bottom:1.1px solid #000!important;font-size:{sz}!important;}}")
EXTRACSS = {
    "小抄": _NOBAR.format(sz="8.6px"),
    "作业讲解": _NOBAR.format(sz="8px"),
}

TEMPLATE = r"""<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<title>__TITLE__</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js"
 onload="renderMathInElement(document.body,{delimiters:[{left:'$$',right:'$$',display:true},{left:'$',right:'$',display:false}],throwOnError:false});document.body.dataset.ready='1';"></script>
<style>
@page { size: A4; margin: __MARGIN__; }
* { box-sizing: border-box; }
html,body { margin:0; padding:0; }
body { font-family:-apple-system,"PingFang SC","Heiti SC","Songti SC",serif;
  font-size:__FONT__; line-height:__LH__; color:#000;
  column-count:__COLS__; column-gap:4mm; column-fill:auto; }
h1 { font-size:11.5px; column-span:all; margin:0 0 2px; padding-bottom:1px;
  border-bottom:1.6px solid #000; text-align:center; }
h2 { font-size:8.8px; margin:4px 0 1px; padding:0.5px 3px; background:#000; color:#fff;
  break-after:avoid; break-inside:avoid; }
h3 { font-size:7.9px; font-weight:700; margin:2.5px 0 0.5px; border-bottom:0.5px solid #888;
  break-after:avoid; }
p { margin:0.6px 0; }
ul,ol { margin:0.6px 0; padding-left:1.1em; }
li { margin:0.3px 0; }
b,strong { font-weight:700; }
blockquote { margin:0.6px 0; padding:0.3px 0 0.3px 3px; border-left:2px solid #aaa;
  color:#1a1a1a; background:#f4f4f4; }
hr { border:0; border-top:0.5px dashed #bbb; margin:2px 0; }
table { border-collapse:collapse; width:100%; font-size:0.9em; margin:1.5px 0;
  break-inside:avoid; table-layout:auto; }
th,td { border:0.4px solid #555; padding:0 1.5px; vertical-align:top; word-break:break-word; }
th { background:#e2e2e2; font-weight:700; }
pre { font-size:5.2px; line-height:1.05; white-space:pre; overflow:hidden;
  margin:1px 0; padding:1px; background:#f4f4f4; border:0.4px solid #ccc; break-inside:avoid; }
code { font-family:Menlo,monospace; }
.katex { font-size:1em; }
.katex-display { margin:1px 0; }
.katex-display>.katex { white-space:normal; }
h2,h3,table,pre,tr { break-inside:avoid; }
__EXTRACSS__
</style></head><body>
__BODY__
</body></html>"""

CODE, MB, MI = {}, {}, {}

def protect(text):
    def repl_code(m):
        k = "\x00C%d\x00" % len(CODE)
        body = m.group(1)
        body = re.sub(r"^[a-zA-Z0-9]*\n", "", body)  # 去掉语言标记行
        CODE[k] = "<pre><code>" + html.escape(body) + "</code></pre>"
        return "\n" + k + "\n"
    text = re.sub(r"```(.*?)```", repl_code, text, flags=re.S)
    def repl_mb(m):
        k = "\x00B%d\x00" % len(MB); MB[k] = m.group(0); return k
    text = re.sub(r"\$\$.*?\$\$", repl_mb, text, flags=re.S)
    def repl_mi(m):
        k = "\x00I%d\x00" % len(MI); MI[k] = m.group(0); return k
    text = re.sub(r"\$[^$\n]*?\$", repl_mi, text)
    return text

def inline(s):
    s = html.escape(s, quote=False)
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"~~(.+?)~~", r"<s>\1</s>", s)
    s = re.sub(r"\[\[(.+?)\]\]", r"\1", s)
    s = re.sub(r"`([^`]+?)`", r"<code>\1</code>", s)
    return s

def split_row(line):
    line = line.strip()
    if line.startswith("|"): line = line[1:]
    if line.endswith("|"): line = line[:-1]
    return [c.strip() for c in line.split("|")]

def is_sep(line):
    return bool(re.match(r"^\s*\|?\s*:?-{2,}.*\|", line)) and set(line.strip()) <= set("|:- ")

def convert(md):
    # 去 YAML frontmatter
    lines = md.split("\n")
    if lines and lines[0].strip() == "---":
        for i in range(1, len(lines)):
            if lines[i].strip() == "---":
                lines = lines[i+1:]; break
    md = protect("\n".join(lines))
    lines = md.split("\n")
    out, i, n = [], 0, len(lines)
    para = []
    def flush():
        if para:
            out.append("<p>" + inline(" ".join(para)) + "</p>"); para.clear()
    while i < n:
        ln = lines[i]
        st = ln.strip()
        if not st:
            flush(); i += 1; continue
        if st.startswith("\x00C"):  # 代码占位
            flush(); out.append(st); i += 1; continue
        m = re.match(r"^(#{1,4})\s+(.*)$", st)
        if m:
            flush(); lvl = len(m.group(1)); out.append("<h%d>%s</h%d>" % (lvl, inline(m.group(2)), lvl)); i += 1; continue
        if st in ("---", "***", "___"):
            flush(); out.append("<hr>"); i += 1; continue
        # 表格
        if st.startswith("|") and i+1 < n and is_sep(lines[i+1]):
            flush()
            head = split_row(lines[i]); i += 2
            rows = []
            while i < n and lines[i].strip().startswith("|"):
                rows.append(split_row(lines[i])); i += 1
            t = ["<table><thead><tr>"] + ["<th>%s</th>" % inline(c) for c in head] + ["</tr></thead><tbody>"]
            for r in rows:
                t.append("<tr>" + "".join("<td>%s</td>" % inline(c) for c in r) + "</tr>")
            t.append("</tbody></table>")
            out.append("".join(t)); continue
        # 引用块
        if st.startswith(">"):
            flush(); buf = []
            while i < n and lines[i].strip().startswith(">"):
                buf.append(inline(re.sub(r"^\s*>+\s?", "", lines[i]))); i += 1
            out.append("<blockquote>" + "<br>".join(buf) + "</blockquote>"); continue
        # 列表
        if re.match(r"^\s*([-*]|\d+\.)\s+", ln):
            flush(); items = []
            while i < n and re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]):
                raw = lines[i]
                indent = len(raw) - len(raw.lstrip(" "))
                txt = re.sub(r"^\s*([-*]|\d+\.)\s+", "", raw)
                ml = (indent // 2) * 1.1
                items.append('<li style="margin-left:%.1fem">%s</li>' % (ml, inline(txt))); i += 1
            out.append("<ul>" + "".join(items) + "</ul>"); continue
        # 普通段落
        para.append(st); i += 1
    flush()
    htmlout = "\n".join(out)
    # 还原占位
    for k, v in MB.items(): htmlout = htmlout.replace(k, v)
    for k, v in MI.items(): htmlout = htmlout.replace(k, v)
    for k, v in CODE.items(): htmlout = htmlout.replace(k, v)
    return htmlout

for name in SRC:
    CODE.clear(); MB.clear(); MI.clear()
    with open(SRC[name], encoding="utf-8") as f:
        md = f.read()
    body = convert(md)
    page = (TEMPLATE.replace("__TITLE__", name)
            .replace("__FONT__", FONT[name]).replace("__COLS__", COLS[name])
            .replace("__LH__", LH[name]).replace("__MARGIN__", MARGIN[name])
            .replace("__EXTRACSS__", EXTRACSS[name])
            .replace("__BODY__", body))
    with open(OUT[name], "w", encoding="utf-8") as f:
        f.write(page)
    print("written:", OUT[name], "chars:", len(page))
