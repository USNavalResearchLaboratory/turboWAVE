" Syntax highlighting for turboWAVE input files:
" Put tw3d.vim into ~/.vim/syntax/ and put filetype.vim into ~/.vim
" This causes stdin, stdin.txt, and *.tw3d to be highlighted
" To highlight any other file from within vim, use ``:set syntax=tw3d``

if exists("did_load_filetypes")
	finish
endif
augroup filetypedetect
	au! BufRead,BufNewFile stdin setfiletype tw3d
	au! BufRead,BufNewFile stdin.txt setfiletype tw3d
	au! BufRead,BufNewFile *.tw3d setfiletype tw3d
augroup END
