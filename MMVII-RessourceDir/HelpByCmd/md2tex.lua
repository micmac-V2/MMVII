-- Pandoc filter used by HelpByCmd.cmake : Markdown help of a command -> body of the MMVII LaTeX help.

-- Drop the title block at the beginning of the file, whatever its number of lines :
-- a "::: center" div, or a level 1 header "Help for command ..." (the title is made by \CmdTitle)
function Pandoc(doc)
  local blocks = doc.blocks
  while #blocks > 0 do
    local b = blocks[1]
    local isCenter = (b.t == "Div") and b.classes:includes("center")
    local isTitle  = (b.t == "Header") and (b.level == 1) and pandoc.utils.stringify(b):match("^Help for command")
    if not (isCenter or isTitle) then break end
    table.remove(blocks, 1)
  end
  doc.blocks = blocks
  return doc
end

-- First section : macro of MMVII.sty
function Header(h)
  if (h.level == 1) and pandoc.utils.stringify(h):match("^Descri") then
    return pandoc.RawBlock("latex", "\\SecDescr")
  end
end

-- Inline code : style of the names of files and parameters
function Code(c)
  local tex = pandoc.write(pandoc.Pandoc({pandoc.Plain({c})}), "latex")
  tex = tex:gsub("^%s*\\texttt", "\\FileSty"):gsub("%s+$", "")
  return pandoc.RawInline("latex", tex)
end
