" Syntax highlighting for turboWAVE input files:
" Put turbowave.vim into ~/.vim/syntax/ and put filetype.vim into ~/.vim
" This causes stdin, stdin.txt, and *.tw to be highlighted
" To highlight any other file from within vim, use ``:set syntax=turbowave``

if exists("did_load_filetypes")
	finish
endif
augroup filetypedetect
	au! BufRead,BufNewFile stdin setfiletype turbowave
	au! BufRead,BufNewFile stdin.txt setfiletype turbowave
	au! BufRead,BufNewFile *.tw setfiletype turbowave
augroup END
