module.exports = grammar(
{
  name: 'turbowave',

  extras: $ => [/\s|,/,$.comment,],

  rules:
  {
    input_file: $ => repeat($._top),

    _top: $ => choice($.include,$.define,$._directive),

    _nested: $ => choice($.include,$.define,$._nested_directive),

    include: $ => seq('#include',choice($.identifier,$.string_literal)),

    define: $ => seq('#define',$.define_key,$.define_value),

    _directive: $ => choice($.assignment,$._statement),

    _nested_directive: $ => choice($.assignment,$._nested_statement),

    assignment: $ => seq($.identifier_sequence,'=',choice($.value,$.tuple,$.list)),

    _statement: $ => choice($.new,$.associative_new,$.generate,$.reaction,$.collision,$.excitation),

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
      seq('for',choice($.string_literal,$.identifier)),
      $.block),

    generate: $ => seq(
      'generate',
      $.identifier_sequence,
      optional($.string_literal),
      $.block),

    reaction: $ => seq('new','reaction','=',$.full_formula,$.rate,$.range),

    collision: $ => seq(
      'new','collision','=',
      $.identifier,'<->',$.identifier,
      choice('coulomb',seq('cross','section','=',$.decimal))),

    excitation: $ => seq(
      'new','excitation','=',
      $.identifier,'->',$.identifier,
      'level','=',$.decimal,
      $.rate),

    get: $ => seq('get',choice($.identifier,$.string_literal)),

    define_key: $ => /\$[a-zA-Z_]\w*/,
    identifier: $ => /[1-3]?[a-zA-Z_][\w\[\]\+\-\^\.]*/,
    identifier_sequence: $ => repeat1($.identifier),
    _string_literal_single: $ => seq('\'',$.identifier,'\''),
    _string_literal_double: $ => seq('\"',$.identifier,'\"'),
    string_literal: $ => choice($._string_literal_double,$._string_literal_single),
    decimal: $ => /(\+|\-)?([0-9]+\.?[0-9]*|\.[0-9]+)([eE](\+|\-)?[0-9]+)?/,
    boolean: $ => choice('true','false','yes','no','on','off'),
    unit: $ => choice('[deg]','[rad]','[mrad]','[urad]','[cm2]','[m2]','[cm2/s]','[m2/s]','[um]','[mm]','[cm]','[m]','[fs]','[ps]','[ns]','[us]','[s]','[/m3]','[/cm3]','[J/m3]','[J/cm3]','[eV]','[K]','[V]','[webers/m]','[G*cm]','[V/m]','[V/cm]','[T]','[G]'),
    dimension: $ => seq($.decimal,$.unit),
    value: $ => choice($.identifier,$.decimal,$.dimension,$.define_key,$.string_literal,$.boolean),
    define_value: $ => choice($.identifier,$.decimal,$.dimension,$.boolean),
    block: $ => seq('{',repeat($._nested),'}'),
    tuple: $ => seq('(',repeat1($.value),')'),
    list: $ => seq('{',repeat1($.value),'}'),

    comment: $ => token(choice(
        seq('//', /(\\(.|\r?\n)|[^\\\n])*/),
        seq('/*',/[^*]*\*+([^/*][^*]*\*+)*/,'/'))),

    // Reaction rule is complex, broken down into the following pieces:
    _pterm: $ => seq(/\s+\+\s+/,choice($.identifier,$.decimal,$.define_key)),
    _nterm: $ => seq(/\s+\-\s+/,choice($.decimal,$.define_key)),
    chems: $ => seq(choice($.identifier,$.define_key),repeat(choice($._pterm,$._nterm))),
    sub_formula: $ => seq($.chems,'->',$.chems),
    full_formula: $ => seq('\{',$.sub_formula,repeat(seq(':',$.sub_formula)),'\}'),
    arrhenius: $ => seq('rate','=',$.decimal,$.decimal,$.decimal),
    janev: $ => seq('janev_rate','=',$.decimal,$.decimal,$.decimal,$.decimal,$.decimal,$.decimal,$.decimal,$.decimal,$.decimal),
    rate: $=> choice($.arrhenius,$.janev),
    range: $=> seq($.identifier,'(',optional($.decimal),':',optional($.decimal),')'),

  }
});
