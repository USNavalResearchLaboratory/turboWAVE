lua <<EOF
local parser_config = require "nvim-treesitter.parsers".get_parser_configs()
parser_config.turbowave = {
	install_info = {
		url = "https://github.com/dfgordon/tree-sitter-turbowave",
		files = {"src/parser.c"}
	},
}
require'nvim-treesitter.configs'.setup {
  ensure_installed = "all",
  sync_install = false,
  ignore_install = { "javascript" },
  highlight = {
  	enable = true,
  	disable = { },
  	additional_vim_regex_highlighting = false,
  },
}
EOF
