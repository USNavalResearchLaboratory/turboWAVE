" Put this in ~/.vim if you want to use java syntax highlights
" when editing the input file stdin using vim editor
if exists("did_load_filetypes")
	finish
endif
augroup filetypedetect
	au! BufRead,BufNewFile stdin	setfiletype java
	au! BufRead,BufNewFile stdin.txt	setfiletype java
augroup END
