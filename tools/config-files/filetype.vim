" Syntax highlighting for turboWAVE input files:
" Put tw3d.vim into ~/.vim/synax/ to highlight files with .tw3d extension
" ALSO put filetype.vim into ~/.vim to highlight stdin and stdin.txt

if exists("did_load_filetypes")
	finish
endif
augroup filetypedetect
	au! BufRead,BufNewFile stdin	setfiletype tw3d
	au! BufRead,BufNewFile stdin.txt	setfiletype tw3d
augroup END
