" Vim syntax file
" Language:	turboWAVE
" Maintainer: Daniel Gordon

" Syntax highlighting for turboWAVE input files:
" Put turbowave.vim into ~/.vim/syntax/ and put filetype.vim into ~/.vim
" This causes stdin, stdin.txt, and *.tw to be highlighted
" To highlight any other file from within vim, use ``:set syntax=turbowave``

" Quit when a (custom) syntax file was already loaded
if exists("b:current_syntax")
  finish
endif

syn match twDefine "#define"
syn match twInclude "#include"
syn keyword	twStatement	new generate

syn match twUnit  "\v(\s|\(|\{|\,|\=)\%[0-9]+\.=[0-9]*[eE]=[\+\-]=[0-9]*(um|mm|cm|m|fs|ps|ns|us|s|m\-3|cm\-3|Jm3|Jcm3|eV|K)(\s|\)|\}|\,|\n)"hs=s+1,he=e-1
syn match twUnit2 "\v(\s|\(|\{|\,|\=)\%\.[0-9]+[eE]=[\+\-][0-9]*(um|mm|cm|m|fs|ps|ns|us|s|m\-3|cm\-3|Jm3|Jcm3|eV|K)(\s|\)|\}|\,|\n)"hs=s+1,he=e-1
syn region twUserMacro display start="\$" end="\s\|$"he=s-1
syn region twComment start="/\*" end="\*/"
syn region twCommentL start="//" end="$" keepend

syn keyword twconst none full
syn keyword twload deterministic statistical variable fixed triggered maintained
syn keyword twbool true false on off yes no
syn keyword twboundary absorbing emitting reflecting axisymmetric periodic
syn keyword twmode airy_disc hermite laguerre bessel plane multipole
syn keyword twioniz adk ppt mpi
syn keyword twshape quintic sin2 sech

syn region twString1 start="'" end="'"
syn region twString2 start='"' end='"'

hi def link twtemp Statement
hi def link twDefine	Define
hi def link twInclude	Define
hi def link twStatement		Statement
hi def link twLabel			Label
hi def link twUnit			Macro
hi def link twUnit2			Macro
hi def link twUserMacro		Macro
hi def link twComment		Comment
hi def link twCommentL		Comment
hi def link twconst		Constant
hi def link twload		Constant
hi def link twbool		Constant
hi def link twboundary		Constant
hi def link twmode		Constant
hi def link twioniz		Constant
hi def link twshape		Constant
hi def link twString1		String
hi def link twString2		String
hi def link twCurlyError	Error

let b:current_syntax = "tw3d"
