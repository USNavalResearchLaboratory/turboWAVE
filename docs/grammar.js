module.exports = grammar(
{
	name: 'turbowave',

	extras: $ => [/\s|\\\r?\n/,',',$.comment,],

	rules:
	{
		input_file: $ => repeat($._top),

		_top: $ => choice($.include,$.define,$._directive),

		_nested: $ => choice($.include,$.define,$._nested_directive),

		include: $ => seq('#include',$.identifier),

		define: $ => seq('#define',$.define_key,$.define_value),

		_directive: $ => choice($.assignment,$._statement),

		_nested_directive: $ => choice($.assignment,$._nested_statement),

		assignment: $ => seq($.identifier_sequence,'=',choice($.value,$.tuple,$.list)),

		_statement: $ => choice($.new,$.associative_new,$.generate),

		_nested_statement: $ => choice($.get,$.new),

		new: $ => seq(
			'new',
			$.identifier_sequence,
			optional($.string_literal),
			$.block),

		associative_new: $ => seq(
			'new',
			$.identifier_sequence,
			optional($.string_literal),
			seq('for',$.string_literal),
			$.block),

		generate: $ => seq(
			'generate',
			$.identifier_sequence,
			$.string_literal,
			$.block),

		get: $ => seq('get',choice($.identifier,$.string_literal)),

		define_key: $ => /\$[a-zA-Z_]\w*/,
		identifier: $ => /[a-zA-Z_]\w*/,
		identifier_sequence: $ => repeat1($.identifier),
		string_literal: $ => /(\'[a-zA-Z_]\w*\'|\"[a-zA-Z_]\w*\")/,
		decimal1: $ => /[0-9]+\.?[0-9]*[eE]?(\+|\-)?[0-9]*/,
		decimal2: $ => /\.[0-9]+[eE]?(\+|\-)?[0-9]*/,

		value: $ => choice($.identifier,$.decimal1,$.decimal2,$.unit_macro,$.define_key),

		define_value: $ => choice($.identifier,$.decimal1,$.decimal2,$.unit_macro),

		unit_macro: $ => seq(optional('-'),'%',choice($.decimal1,$.decimal2),$.unit),

		unit: $ => choice('deg','rad','mrad','urad','cm2','m2','cm2s','m2s','um','mm','cm','m','fs','ps','ns','us','s','m-3','cm-3','Jm3','Jcm3','eV','K'),

		block: $ => seq('{',repeat($._nested),'}'),

		tuple: $ => seq('(',repeat1($.value),')'),

		list: $ => seq('{',repeat1($.value),'}'),

		comment: $ => token(choice(
  			seq('//', /(\\(.|\r?\n)|[^\\\n])*/),
  			seq('/*',/[^*]*\*+([^/*][^*]*\*+)*/,'/'))),
	}
});
