#include <parser.h>

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif

#ifdef _MSC_VER
#pragma optimize("", off)
#elif defined(__clang__)
#pragma clang optimize off
#elif defined(__GNUC__)
#pragma GCC optimize ("O0")
#endif

#define LANGUAGE_VERSION 13
#define STATE_COUNT 338
#define LARGE_STATE_COUNT 2
#define SYMBOL_COUNT 137
#define ALIAS_COUNT 0
#define TOKEN_COUNT 88
#define EXTERNAL_TOKEN_COUNT 0
#define FIELD_COUNT 0
#define MAX_ALIAS_SEQUENCE_LENGTH 16
#define PRODUCTION_ID_COUNT 1

enum {
  anon_sym_POUNDinclude = 1,
  anon_sym_POUNDdefine = 2,
  anon_sym_POUNDifdef = 3,
  anon_sym_POUNDifndef = 4,
  anon_sym_POUNDendif = 5,
  anon_sym_EQ = 6,
  anon_sym_new = 7,
  anon_sym_for = 8,
  anon_sym_generate = 9,
  anon_sym_reaction = 10,
  anon_sym_collision = 11,
  anon_sym_LT_DASH_GT = 12,
  anon_sym_coulomb = 13,
  anon_sym_cross = 14,
  anon_sym_section = 15,
  anon_sym_metallic = 16,
  anon_sym_ks = 17,
  anon_sym_fermi_energy_ev = 18,
  anon_sym_ref_density = 19,
  anon_sym_excitation = 20,
  anon_sym_DASH_GT = 21,
  anon_sym_level = 22,
  anon_sym_get = 23,
  anon_sym_POUNDelse = 24,
  sym_define_key = 25,
  sym_define_ref = 26,
  sym_decimal = 27,
  anon_sym_true = 28,
  anon_sym_false = 29,
  anon_sym_yes = 30,
  anon_sym_no = 31,
  anon_sym_on = 32,
  anon_sym_off = 33,
  anon_sym_LBRACKdeg_RBRACK = 34,
  anon_sym_LBRACKrad_RBRACK = 35,
  anon_sym_LBRACKmrad_RBRACK = 36,
  anon_sym_LBRACKurad_RBRACK = 37,
  anon_sym_LBRACKum_RBRACK = 38,
  anon_sym_LBRACKmm_RBRACK = 39,
  anon_sym_LBRACKcm_RBRACK = 40,
  anon_sym_LBRACKm_RBRACK = 41,
  anon_sym_LBRACKfs_RBRACK = 42,
  anon_sym_LBRACKps_RBRACK = 43,
  anon_sym_LBRACKns_RBRACK = 44,
  anon_sym_LBRACKus_RBRACK = 45,
  anon_sym_LBRACKs_RBRACK = 46,
  anon_sym_LBRACK_SLASHm3_RBRACK = 47,
  anon_sym_LBRACK_SLASHcm3_RBRACK = 48,
  anon_sym_LBRACKkg_SLASHm3_RBRACK = 49,
  anon_sym_LBRACKg_SLASHcm3_RBRACK = 50,
  anon_sym_LBRACKJ_SLASHm3_RBRACK = 51,
  anon_sym_LBRACKJ_SLASHcm3_RBRACK = 52,
  anon_sym_LBRACKeV_RBRACK = 53,
  anon_sym_LBRACKK_RBRACK = 54,
  anon_sym_LBRACKPa_RBRACK = 55,
  anon_sym_LBRACKdynes_SLASHcm2_RBRACK = 56,
  anon_sym_LBRACKbar_RBRACK = 57,
  anon_sym_LBRACKergs_SLASHg_RBRACK = 58,
  anon_sym_LBRACKJ_SLASHkg_RBRACK = 59,
  anon_sym_LBRACKcm2_RBRACK = 60,
  anon_sym_LBRACKm2_RBRACK = 61,
  anon_sym_LBRACKcm2_SLASHs_RBRACK = 62,
  anon_sym_LBRACKm2_SLASHs_RBRACK = 63,
  anon_sym_LBRACKV_RBRACK = 64,
  anon_sym_LBRACKwebers_SLASHm_RBRACK = 65,
  anon_sym_LBRACKG_STARcm_RBRACK = 66,
  anon_sym_LBRACKV_SLASHm_RBRACK = 67,
  anon_sym_LBRACKV_SLASHcm_RBRACK = 68,
  anon_sym_LBRACKT_RBRACK = 69,
  anon_sym_LBRACKG_RBRACK = 70,
  sym_identifier = 71,
  anon_sym_1d = 72,
  anon_sym_2d = 73,
  anon_sym_3d = 74,
  anon_sym_SQUOTE = 75,
  anon_sym_DQUOTE = 76,
  anon_sym_LBRACE = 77,
  anon_sym_RBRACE = 78,
  anon_sym_LPAREN = 79,
  anon_sym_RPAREN = 80,
  aux_sym__define_value_token1 = 81,
  sym_comment = 82,
  aux_sym__pterm_token1 = 83,
  aux_sym__nterm_token1 = 84,
  anon_sym_COLON = 85,
  anon_sym_rate = 86,
  anon_sym_janev_rate = 87,
  sym_input_file = 88,
  sym__top = 89,
  sym__nested = 90,
  sym_include = 91,
  sym_define = 92,
  sym_ifxdef = 93,
  sym__directive = 94,
  sym__nested_directive = 95,
  sym_assignment = 96,
  sym__statement = 97,
  sym__nested_statement = 98,
  sym_new = 99,
  sym_associative_new = 100,
  sym_generate = 101,
  sym_reaction = 102,
  sym_collision = 103,
  sym_excitation = 104,
  sym_get = 105,
  sym_else_block = 106,
  sym_boolean = 107,
  sym_unit = 108,
  sym_dimension = 109,
  sym_special_keys = 110,
  sym_obj_key = 111,
  sym__string_literal_single = 112,
  sym__string_literal_double = 113,
  sym_string_literal = 114,
  sym__value = 115,
  sym__rawqty = 116,
  sym_block = 117,
  sym_tuple = 118,
  sym_list = 119,
  sym__define_value = 120,
  sym__pterm = 121,
  sym__nterm = 122,
  sym_chems = 123,
  sym_sub_formula = 124,
  sym_full_formula = 125,
  sym_arrhenius = 126,
  sym_janev = 127,
  sym_rate = 128,
  sym_range = 129,
  aux_sym_input_file_repeat1 = 130,
  aux_sym_obj_key_repeat1 = 131,
  aux_sym_block_repeat1 = 132,
  aux_sym_tuple_repeat1 = 133,
  aux_sym__define_value_repeat1 = 134,
  aux_sym_chems_repeat1 = 135,
  aux_sym_full_formula_repeat1 = 136,
};

static const char * const ts_symbol_names[] = {
  [ts_builtin_sym_end] = "end",
  [anon_sym_POUNDinclude] = "#include",
  [anon_sym_POUNDdefine] = "#define",
  [anon_sym_POUNDifdef] = "#ifdef",
  [anon_sym_POUNDifndef] = "#ifndef",
  [anon_sym_POUNDendif] = "#endif",
  [anon_sym_EQ] = "=",
  [anon_sym_new] = "new",
  [anon_sym_for] = "for",
  [anon_sym_generate] = "generate",
  [anon_sym_reaction] = "reaction",
  [anon_sym_collision] = "collision",
  [anon_sym_LT_DASH_GT] = "<->",
  [anon_sym_coulomb] = "coulomb",
  [anon_sym_cross] = "cross",
  [anon_sym_section] = "section",
  [anon_sym_metallic] = "metallic",
  [anon_sym_ks] = "ks",
  [anon_sym_fermi_energy_ev] = "fermi_energy_ev",
  [anon_sym_ref_density] = "ref_density",
  [anon_sym_excitation] = "excitation",
  [anon_sym_DASH_GT] = "->",
  [anon_sym_level] = "level",
  [anon_sym_get] = "get",
  [anon_sym_POUNDelse] = "#else",
  [sym_define_key] = "define_key",
  [sym_define_ref] = "define_ref",
  [sym_decimal] = "decimal",
  [anon_sym_true] = "true",
  [anon_sym_false] = "false",
  [anon_sym_yes] = "yes",
  [anon_sym_no] = "no",
  [anon_sym_on] = "on",
  [anon_sym_off] = "off",
  [anon_sym_LBRACKdeg_RBRACK] = "[deg]",
  [anon_sym_LBRACKrad_RBRACK] = "[rad]",
  [anon_sym_LBRACKmrad_RBRACK] = "[mrad]",
  [anon_sym_LBRACKurad_RBRACK] = "[urad]",
  [anon_sym_LBRACKum_RBRACK] = "[um]",
  [anon_sym_LBRACKmm_RBRACK] = "[mm]",
  [anon_sym_LBRACKcm_RBRACK] = "[cm]",
  [anon_sym_LBRACKm_RBRACK] = "[m]",
  [anon_sym_LBRACKfs_RBRACK] = "[fs]",
  [anon_sym_LBRACKps_RBRACK] = "[ps]",
  [anon_sym_LBRACKns_RBRACK] = "[ns]",
  [anon_sym_LBRACKus_RBRACK] = "[us]",
  [anon_sym_LBRACKs_RBRACK] = "[s]",
  [anon_sym_LBRACK_SLASHm3_RBRACK] = "[/m3]",
  [anon_sym_LBRACK_SLASHcm3_RBRACK] = "[/cm3]",
  [anon_sym_LBRACKkg_SLASHm3_RBRACK] = "[kg/m3]",
  [anon_sym_LBRACKg_SLASHcm3_RBRACK] = "[g/cm3]",
  [anon_sym_LBRACKJ_SLASHm3_RBRACK] = "[J/m3]",
  [anon_sym_LBRACKJ_SLASHcm3_RBRACK] = "[J/cm3]",
  [anon_sym_LBRACKeV_RBRACK] = "[eV]",
  [anon_sym_LBRACKK_RBRACK] = "[K]",
  [anon_sym_LBRACKPa_RBRACK] = "[Pa]",
  [anon_sym_LBRACKdynes_SLASHcm2_RBRACK] = "[dynes/cm2]",
  [anon_sym_LBRACKbar_RBRACK] = "[bar]",
  [anon_sym_LBRACKergs_SLASHg_RBRACK] = "[ergs/g]",
  [anon_sym_LBRACKJ_SLASHkg_RBRACK] = "[J/kg]",
  [anon_sym_LBRACKcm2_RBRACK] = "[cm2]",
  [anon_sym_LBRACKm2_RBRACK] = "[m2]",
  [anon_sym_LBRACKcm2_SLASHs_RBRACK] = "[cm2/s]",
  [anon_sym_LBRACKm2_SLASHs_RBRACK] = "[m2/s]",
  [anon_sym_LBRACKV_RBRACK] = "[V]",
  [anon_sym_LBRACKwebers_SLASHm_RBRACK] = "[webers/m]",
  [anon_sym_LBRACKG_STARcm_RBRACK] = "[G*cm]",
  [anon_sym_LBRACKV_SLASHm_RBRACK] = "[V/m]",
  [anon_sym_LBRACKV_SLASHcm_RBRACK] = "[V/cm]",
  [anon_sym_LBRACKT_RBRACK] = "[T]",
  [anon_sym_LBRACKG_RBRACK] = "[G]",
  [sym_identifier] = "identifier",
  [anon_sym_1d] = "1d",
  [anon_sym_2d] = "2d",
  [anon_sym_3d] = "3d",
  [anon_sym_SQUOTE] = "'",
  [anon_sym_DQUOTE] = "\"",
  [anon_sym_LBRACE] = "{",
  [anon_sym_RBRACE] = "}",
  [anon_sym_LPAREN] = "(",
  [anon_sym_RPAREN] = ")",
  [aux_sym__define_value_token1] = "_define_value_token1",
  [sym_comment] = "comment",
  [aux_sym__pterm_token1] = "_pterm_token1",
  [aux_sym__nterm_token1] = "_nterm_token1",
  [anon_sym_COLON] = ":",
  [anon_sym_rate] = "rate",
  [anon_sym_janev_rate] = "janev_rate",
  [sym_input_file] = "input_file",
  [sym__top] = "_top",
  [sym__nested] = "_nested",
  [sym_include] = "include",
  [sym_define] = "define",
  [sym_ifxdef] = "ifxdef",
  [sym__directive] = "_directive",
  [sym__nested_directive] = "_nested_directive",
  [sym_assignment] = "assignment",
  [sym__statement] = "_statement",
  [sym__nested_statement] = "_nested_statement",
  [sym_new] = "new",
  [sym_associative_new] = "associative_new",
  [sym_generate] = "generate",
  [sym_reaction] = "reaction",
  [sym_collision] = "collision",
  [sym_excitation] = "excitation",
  [sym_get] = "get",
  [sym_else_block] = "else_block",
  [sym_boolean] = "boolean",
  [sym_unit] = "unit",
  [sym_dimension] = "dimension",
  [sym_special_keys] = "special_keys",
  [sym_obj_key] = "obj_key",
  [sym__string_literal_single] = "_string_literal_single",
  [sym__string_literal_double] = "_string_literal_double",
  [sym_string_literal] = "string_literal",
  [sym__value] = "_value",
  [sym__rawqty] = "_rawqty",
  [sym_block] = "block",
  [sym_tuple] = "tuple",
  [sym_list] = "list",
  [sym__define_value] = "_define_value",
  [sym__pterm] = "_pterm",
  [sym__nterm] = "_nterm",
  [sym_chems] = "chems",
  [sym_sub_formula] = "sub_formula",
  [sym_full_formula] = "full_formula",
  [sym_arrhenius] = "arrhenius",
  [sym_janev] = "janev",
  [sym_rate] = "rate",
  [sym_range] = "range",
  [aux_sym_input_file_repeat1] = "input_file_repeat1",
  [aux_sym_obj_key_repeat1] = "obj_key_repeat1",
  [aux_sym_block_repeat1] = "block_repeat1",
  [aux_sym_tuple_repeat1] = "tuple_repeat1",
  [aux_sym__define_value_repeat1] = "_define_value_repeat1",
  [aux_sym_chems_repeat1] = "chems_repeat1",
  [aux_sym_full_formula_repeat1] = "full_formula_repeat1",
};

static const TSSymbol ts_symbol_map[] = {
  [ts_builtin_sym_end] = ts_builtin_sym_end,
  [anon_sym_POUNDinclude] = anon_sym_POUNDinclude,
  [anon_sym_POUNDdefine] = anon_sym_POUNDdefine,
  [anon_sym_POUNDifdef] = anon_sym_POUNDifdef,
  [anon_sym_POUNDifndef] = anon_sym_POUNDifndef,
  [anon_sym_POUNDendif] = anon_sym_POUNDendif,
  [anon_sym_EQ] = anon_sym_EQ,
  [anon_sym_new] = anon_sym_new,
  [anon_sym_for] = anon_sym_for,
  [anon_sym_generate] = anon_sym_generate,
  [anon_sym_reaction] = anon_sym_reaction,
  [anon_sym_collision] = anon_sym_collision,
  [anon_sym_LT_DASH_GT] = anon_sym_LT_DASH_GT,
  [anon_sym_coulomb] = anon_sym_coulomb,
  [anon_sym_cross] = anon_sym_cross,
  [anon_sym_section] = anon_sym_section,
  [anon_sym_metallic] = anon_sym_metallic,
  [anon_sym_ks] = anon_sym_ks,
  [anon_sym_fermi_energy_ev] = anon_sym_fermi_energy_ev,
  [anon_sym_ref_density] = anon_sym_ref_density,
  [anon_sym_excitation] = anon_sym_excitation,
  [anon_sym_DASH_GT] = anon_sym_DASH_GT,
  [anon_sym_level] = anon_sym_level,
  [anon_sym_get] = anon_sym_get,
  [anon_sym_POUNDelse] = anon_sym_POUNDelse,
  [sym_define_key] = sym_define_key,
  [sym_define_ref] = sym_define_ref,
  [sym_decimal] = sym_decimal,
  [anon_sym_true] = anon_sym_true,
  [anon_sym_false] = anon_sym_false,
  [anon_sym_yes] = anon_sym_yes,
  [anon_sym_no] = anon_sym_no,
  [anon_sym_on] = anon_sym_on,
  [anon_sym_off] = anon_sym_off,
  [anon_sym_LBRACKdeg_RBRACK] = anon_sym_LBRACKdeg_RBRACK,
  [anon_sym_LBRACKrad_RBRACK] = anon_sym_LBRACKrad_RBRACK,
  [anon_sym_LBRACKmrad_RBRACK] = anon_sym_LBRACKmrad_RBRACK,
  [anon_sym_LBRACKurad_RBRACK] = anon_sym_LBRACKurad_RBRACK,
  [anon_sym_LBRACKum_RBRACK] = anon_sym_LBRACKum_RBRACK,
  [anon_sym_LBRACKmm_RBRACK] = anon_sym_LBRACKmm_RBRACK,
  [anon_sym_LBRACKcm_RBRACK] = anon_sym_LBRACKcm_RBRACK,
  [anon_sym_LBRACKm_RBRACK] = anon_sym_LBRACKm_RBRACK,
  [anon_sym_LBRACKfs_RBRACK] = anon_sym_LBRACKfs_RBRACK,
  [anon_sym_LBRACKps_RBRACK] = anon_sym_LBRACKps_RBRACK,
  [anon_sym_LBRACKns_RBRACK] = anon_sym_LBRACKns_RBRACK,
  [anon_sym_LBRACKus_RBRACK] = anon_sym_LBRACKus_RBRACK,
  [anon_sym_LBRACKs_RBRACK] = anon_sym_LBRACKs_RBRACK,
  [anon_sym_LBRACK_SLASHm3_RBRACK] = anon_sym_LBRACK_SLASHm3_RBRACK,
  [anon_sym_LBRACK_SLASHcm3_RBRACK] = anon_sym_LBRACK_SLASHcm3_RBRACK,
  [anon_sym_LBRACKkg_SLASHm3_RBRACK] = anon_sym_LBRACKkg_SLASHm3_RBRACK,
  [anon_sym_LBRACKg_SLASHcm3_RBRACK] = anon_sym_LBRACKg_SLASHcm3_RBRACK,
  [anon_sym_LBRACKJ_SLASHm3_RBRACK] = anon_sym_LBRACKJ_SLASHm3_RBRACK,
  [anon_sym_LBRACKJ_SLASHcm3_RBRACK] = anon_sym_LBRACKJ_SLASHcm3_RBRACK,
  [anon_sym_LBRACKeV_RBRACK] = anon_sym_LBRACKeV_RBRACK,
  [anon_sym_LBRACKK_RBRACK] = anon_sym_LBRACKK_RBRACK,
  [anon_sym_LBRACKPa_RBRACK] = anon_sym_LBRACKPa_RBRACK,
  [anon_sym_LBRACKdynes_SLASHcm2_RBRACK] = anon_sym_LBRACKdynes_SLASHcm2_RBRACK,
  [anon_sym_LBRACKbar_RBRACK] = anon_sym_LBRACKbar_RBRACK,
  [anon_sym_LBRACKergs_SLASHg_RBRACK] = anon_sym_LBRACKergs_SLASHg_RBRACK,
  [anon_sym_LBRACKJ_SLASHkg_RBRACK] = anon_sym_LBRACKJ_SLASHkg_RBRACK,
  [anon_sym_LBRACKcm2_RBRACK] = anon_sym_LBRACKcm2_RBRACK,
  [anon_sym_LBRACKm2_RBRACK] = anon_sym_LBRACKm2_RBRACK,
  [anon_sym_LBRACKcm2_SLASHs_RBRACK] = anon_sym_LBRACKcm2_SLASHs_RBRACK,
  [anon_sym_LBRACKm2_SLASHs_RBRACK] = anon_sym_LBRACKm2_SLASHs_RBRACK,
  [anon_sym_LBRACKV_RBRACK] = anon_sym_LBRACKV_RBRACK,
  [anon_sym_LBRACKwebers_SLASHm_RBRACK] = anon_sym_LBRACKwebers_SLASHm_RBRACK,
  [anon_sym_LBRACKG_STARcm_RBRACK] = anon_sym_LBRACKG_STARcm_RBRACK,
  [anon_sym_LBRACKV_SLASHm_RBRACK] = anon_sym_LBRACKV_SLASHm_RBRACK,
  [anon_sym_LBRACKV_SLASHcm_RBRACK] = anon_sym_LBRACKV_SLASHcm_RBRACK,
  [anon_sym_LBRACKT_RBRACK] = anon_sym_LBRACKT_RBRACK,
  [anon_sym_LBRACKG_RBRACK] = anon_sym_LBRACKG_RBRACK,
  [sym_identifier] = sym_identifier,
  [anon_sym_1d] = anon_sym_1d,
  [anon_sym_2d] = anon_sym_2d,
  [anon_sym_3d] = anon_sym_3d,
  [anon_sym_SQUOTE] = anon_sym_SQUOTE,
  [anon_sym_DQUOTE] = anon_sym_DQUOTE,
  [anon_sym_LBRACE] = anon_sym_LBRACE,
  [anon_sym_RBRACE] = anon_sym_RBRACE,
  [anon_sym_LPAREN] = anon_sym_LPAREN,
  [anon_sym_RPAREN] = anon_sym_RPAREN,
  [aux_sym__define_value_token1] = aux_sym__define_value_token1,
  [sym_comment] = sym_comment,
  [aux_sym__pterm_token1] = aux_sym__pterm_token1,
  [aux_sym__nterm_token1] = aux_sym__nterm_token1,
  [anon_sym_COLON] = anon_sym_COLON,
  [anon_sym_rate] = anon_sym_rate,
  [anon_sym_janev_rate] = anon_sym_janev_rate,
  [sym_input_file] = sym_input_file,
  [sym__top] = sym__top,
  [sym__nested] = sym__nested,
  [sym_include] = sym_include,
  [sym_define] = sym_define,
  [sym_ifxdef] = sym_ifxdef,
  [sym__directive] = sym__directive,
  [sym__nested_directive] = sym__nested_directive,
  [sym_assignment] = sym_assignment,
  [sym__statement] = sym__statement,
  [sym__nested_statement] = sym__nested_statement,
  [sym_new] = sym_new,
  [sym_associative_new] = sym_associative_new,
  [sym_generate] = sym_generate,
  [sym_reaction] = sym_reaction,
  [sym_collision] = sym_collision,
  [sym_excitation] = sym_excitation,
  [sym_get] = sym_get,
  [sym_else_block] = sym_else_block,
  [sym_boolean] = sym_boolean,
  [sym_unit] = sym_unit,
  [sym_dimension] = sym_dimension,
  [sym_special_keys] = sym_special_keys,
  [sym_obj_key] = sym_obj_key,
  [sym__string_literal_single] = sym__string_literal_single,
  [sym__string_literal_double] = sym__string_literal_double,
  [sym_string_literal] = sym_string_literal,
  [sym__value] = sym__value,
  [sym__rawqty] = sym__rawqty,
  [sym_block] = sym_block,
  [sym_tuple] = sym_tuple,
  [sym_list] = sym_list,
  [sym__define_value] = sym__define_value,
  [sym__pterm] = sym__pterm,
  [sym__nterm] = sym__nterm,
  [sym_chems] = sym_chems,
  [sym_sub_formula] = sym_sub_formula,
  [sym_full_formula] = sym_full_formula,
  [sym_arrhenius] = sym_arrhenius,
  [sym_janev] = sym_janev,
  [sym_rate] = sym_rate,
  [sym_range] = sym_range,
  [aux_sym_input_file_repeat1] = aux_sym_input_file_repeat1,
  [aux_sym_obj_key_repeat1] = aux_sym_obj_key_repeat1,
  [aux_sym_block_repeat1] = aux_sym_block_repeat1,
  [aux_sym_tuple_repeat1] = aux_sym_tuple_repeat1,
  [aux_sym__define_value_repeat1] = aux_sym__define_value_repeat1,
  [aux_sym_chems_repeat1] = aux_sym_chems_repeat1,
  [aux_sym_full_formula_repeat1] = aux_sym_full_formula_repeat1,
};

static const TSSymbolMetadata ts_symbol_metadata[] = {
  [ts_builtin_sym_end] = {
    .visible = false,
    .named = true,
  },
  [anon_sym_POUNDinclude] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_POUNDdefine] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_POUNDifdef] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_POUNDifndef] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_POUNDendif] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_EQ] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_new] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_for] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_generate] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_reaction] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_collision] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LT_DASH_GT] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_coulomb] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_cross] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_section] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_metallic] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_ks] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_fermi_energy_ev] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_ref_density] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_excitation] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_DASH_GT] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_level] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_get] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_POUNDelse] = {
    .visible = true,
    .named = false,
  },
  [sym_define_key] = {
    .visible = true,
    .named = true,
  },
  [sym_define_ref] = {
    .visible = true,
    .named = true,
  },
  [sym_decimal] = {
    .visible = true,
    .named = true,
  },
  [anon_sym_true] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_false] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_yes] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_no] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_on] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_off] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKdeg_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKrad_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKmrad_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKurad_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKum_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKmm_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKcm_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKm_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKfs_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKps_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKns_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKus_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKs_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACK_SLASHm3_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACK_SLASHcm3_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKkg_SLASHm3_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKg_SLASHcm3_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKJ_SLASHm3_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKJ_SLASHcm3_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKeV_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKK_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKPa_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKdynes_SLASHcm2_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKbar_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKergs_SLASHg_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKJ_SLASHkg_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKcm2_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKm2_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKcm2_SLASHs_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKm2_SLASHs_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKV_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKwebers_SLASHm_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKG_STARcm_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKV_SLASHm_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKV_SLASHcm_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKT_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACKG_RBRACK] = {
    .visible = true,
    .named = false,
  },
  [sym_identifier] = {
    .visible = true,
    .named = true,
  },
  [anon_sym_1d] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_2d] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_3d] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_SQUOTE] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_DQUOTE] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LBRACE] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_RBRACE] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_LPAREN] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_RPAREN] = {
    .visible = true,
    .named = false,
  },
  [aux_sym__define_value_token1] = {
    .visible = false,
    .named = false,
  },
  [sym_comment] = {
    .visible = true,
    .named = true,
  },
  [aux_sym__pterm_token1] = {
    .visible = false,
    .named = false,
  },
  [aux_sym__nterm_token1] = {
    .visible = false,
    .named = false,
  },
  [anon_sym_COLON] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_rate] = {
    .visible = true,
    .named = false,
  },
  [anon_sym_janev_rate] = {
    .visible = true,
    .named = false,
  },
  [sym_input_file] = {
    .visible = true,
    .named = true,
  },
  [sym__top] = {
    .visible = false,
    .named = true,
  },
  [sym__nested] = {
    .visible = false,
    .named = true,
  },
  [sym_include] = {
    .visible = true,
    .named = true,
  },
  [sym_define] = {
    .visible = true,
    .named = true,
  },
  [sym_ifxdef] = {
    .visible = true,
    .named = true,
  },
  [sym__directive] = {
    .visible = false,
    .named = true,
  },
  [sym__nested_directive] = {
    .visible = false,
    .named = true,
  },
  [sym_assignment] = {
    .visible = true,
    .named = true,
  },
  [sym__statement] = {
    .visible = false,
    .named = true,
  },
  [sym__nested_statement] = {
    .visible = false,
    .named = true,
  },
  [sym_new] = {
    .visible = true,
    .named = true,
  },
  [sym_associative_new] = {
    .visible = true,
    .named = true,
  },
  [sym_generate] = {
    .visible = true,
    .named = true,
  },
  [sym_reaction] = {
    .visible = true,
    .named = true,
  },
  [sym_collision] = {
    .visible = true,
    .named = true,
  },
  [sym_excitation] = {
    .visible = true,
    .named = true,
  },
  [sym_get] = {
    .visible = true,
    .named = true,
  },
  [sym_else_block] = {
    .visible = true,
    .named = true,
  },
  [sym_boolean] = {
    .visible = true,
    .named = true,
  },
  [sym_unit] = {
    .visible = true,
    .named = true,
  },
  [sym_dimension] = {
    .visible = true,
    .named = true,
  },
  [sym_special_keys] = {
    .visible = true,
    .named = true,
  },
  [sym_obj_key] = {
    .visible = true,
    .named = true,
  },
  [sym__string_literal_single] = {
    .visible = false,
    .named = true,
  },
  [sym__string_literal_double] = {
    .visible = false,
    .named = true,
  },
  [sym_string_literal] = {
    .visible = true,
    .named = true,
  },
  [sym__value] = {
    .visible = false,
    .named = true,
  },
  [sym__rawqty] = {
    .visible = false,
    .named = true,
  },
  [sym_block] = {
    .visible = true,
    .named = true,
  },
  [sym_tuple] = {
    .visible = true,
    .named = true,
  },
  [sym_list] = {
    .visible = true,
    .named = true,
  },
  [sym__define_value] = {
    .visible = false,
    .named = true,
  },
  [sym__pterm] = {
    .visible = false,
    .named = true,
  },
  [sym__nterm] = {
    .visible = false,
    .named = true,
  },
  [sym_chems] = {
    .visible = true,
    .named = true,
  },
  [sym_sub_formula] = {
    .visible = true,
    .named = true,
  },
  [sym_full_formula] = {
    .visible = true,
    .named = true,
  },
  [sym_arrhenius] = {
    .visible = true,
    .named = true,
  },
  [sym_janev] = {
    .visible = true,
    .named = true,
  },
  [sym_rate] = {
    .visible = true,
    .named = true,
  },
  [sym_range] = {
    .visible = true,
    .named = true,
  },
  [aux_sym_input_file_repeat1] = {
    .visible = false,
    .named = false,
  },
  [aux_sym_obj_key_repeat1] = {
    .visible = false,
    .named = false,
  },
  [aux_sym_block_repeat1] = {
    .visible = false,
    .named = false,
  },
  [aux_sym_tuple_repeat1] = {
    .visible = false,
    .named = false,
  },
  [aux_sym__define_value_repeat1] = {
    .visible = false,
    .named = false,
  },
  [aux_sym_chems_repeat1] = {
    .visible = false,
    .named = false,
  },
  [aux_sym_full_formula_repeat1] = {
    .visible = false,
    .named = false,
  },
};

static const TSSymbol ts_alias_sequences[PRODUCTION_ID_COUNT][MAX_ALIAS_SEQUENCE_LENGTH] = {
  [0] = {0},
};

static const uint16_t ts_non_terminal_alias_map[] = {
  0,
};

static inline bool sym_define_key_character_set_1(int32_t c) {
  return (c < ','
    ? (c < '\r'
      ? (c < '\t'
        ? c == 0
        : c <= '\n')
      : (c <= '\r' || (c < '('
        ? c == ' '
        : c <= ')')))
    : (c <= ',' || (c < '{'
      ? (c < '='
        ? c == ':'
        : c <= '=')
      : (c <= '{' || c == '}'))));
}

static bool ts_lex(TSLexer *lexer, TSStateId state) {
  START_LEXER();
  eof = lexer->eof(lexer);
  switch (state) {
    case 0:
      if (eof) ADVANCE(255);
      if (lookahead == '"') ADVANCE(399);
      if (lookahead == '#') ADVANCE(106);
      if (lookahead == '$') ADVANCE(251);
      if (lookahead == '\'') ADVANCE(398);
      if (lookahead == '(') ADVANCE(402);
      if (lookahead == ')') ADVANCE(403);
      if (lookahead == '+') ADVANCE(13);
      if (lookahead == ',') SKIP(0)
      if (lookahead == '-') ADVANCE(11);
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(289);
      if (lookahead == '2') ADVANCE(290);
      if (lookahead == '3') ADVANCE(291);
      if (lookahead == ':') ADVANCE(410);
      if (lookahead == '<') ADVANCE(19);
      if (lookahead == '=') ADVANCE(261);
      if (lookahead == '[') ADVANCE(20);
      if (lookahead == 'c') ADVANCE(198);
      if (lookahead == 'e') ADVANCE(244);
      if (lookahead == 'f') ADVANCE(76);
      if (lookahead == 'g') ADVANCE(115);
      if (lookahead == 'j') ADVANCE(77);
      if (lookahead == 'k') ADVANCE(213);
      if (lookahead == 'l') ADVANCE(116);
      if (lookahead == 'm') ADVANCE(117);
      if (lookahead == 'n') ADVANCE(118);
      if (lookahead == 'o') ADVANCE(143);
      if (lookahead == 'r') ADVANCE(78);
      if (lookahead == 's') ADVANCE(129);
      if (lookahead == 't') ADVANCE(205);
      if (lookahead == 'y') ADVANCE(133);
      if (lookahead == '{') ADVANCE(400);
      if (lookahead == '}') ADVANCE(401);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ') ADVANCE(3);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 1:
      if (lookahead == '\n') ADVANCE(404);
      if (lookahead == '\r') ADVANCE(1);
      if (lookahead == '"') ADVANCE(399);
      if (lookahead == '$') ADVANCE(252);
      if (lookahead == '\'') ADVANCE(398);
      if (lookahead == '(') ADVANCE(402);
      if (lookahead == '+' ||
          lookahead == '-') ADVANCE(13);
      if (lookahead == '\t' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(1)
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(289);
      if (lookahead == '2') ADVANCE(290);
      if (lookahead == '3') ADVANCE(291);
      if (lookahead == '=') ADVANCE(261);
      if (lookahead == '[') ADVANCE(20);
      if (lookahead == 'f') ADVANCE(345);
      if (lookahead == 'g') ADVANCE(350);
      if (lookahead == 'n') ADVANCE(351);
      if (lookahead == 'o') ADVANCE(360);
      if (lookahead == 't') ADVANCE(381);
      if (lookahead == 'y') ADVANCE(355);
      if (lookahead == '{') ADVANCE(400);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 2:
      if (lookahead == '"') ADVANCE(399);
      if (lookahead == '#') ADVANCE(106);
      if (lookahead == '$') ADVANCE(251);
      if (lookahead == '\'') ADVANCE(398);
      if (lookahead == '(') ADVANCE(402);
      if (lookahead == ')') ADVANCE(403);
      if (lookahead == '+') ADVANCE(13);
      if (lookahead == ',') SKIP(2)
      if (lookahead == '-') ADVANCE(11);
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(289);
      if (lookahead == '2') ADVANCE(290);
      if (lookahead == '3') ADVANCE(291);
      if (lookahead == ':') ADVANCE(410);
      if (lookahead == '<') ADVANCE(19);
      if (lookahead == '=') ADVANCE(261);
      if (lookahead == '[') ADVANCE(20);
      if (lookahead == 'c') ADVANCE(198);
      if (lookahead == 'e') ADVANCE(244);
      if (lookahead == 'f') ADVANCE(76);
      if (lookahead == 'g') ADVANCE(115);
      if (lookahead == 'j') ADVANCE(77);
      if (lookahead == 'k') ADVANCE(213);
      if (lookahead == 'l') ADVANCE(116);
      if (lookahead == 'm') ADVANCE(117);
      if (lookahead == 'n') ADVANCE(118);
      if (lookahead == 'o') ADVANCE(143);
      if (lookahead == 'r') ADVANCE(78);
      if (lookahead == 's') ADVANCE(129);
      if (lookahead == 't') ADVANCE(205);
      if (lookahead == 'y') ADVANCE(133);
      if (lookahead == '{') ADVANCE(400);
      if (lookahead == '}') ADVANCE(401);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ') ADVANCE(3);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 3:
      if (lookahead == '"') ADVANCE(399);
      if (lookahead == '#') ADVANCE(106);
      if (lookahead == '$') ADVANCE(251);
      if (lookahead == '\'') ADVANCE(398);
      if (lookahead == '(') ADVANCE(402);
      if (lookahead == ')') ADVANCE(403);
      if (lookahead == '+') ADVANCE(12);
      if (lookahead == ',') SKIP(2)
      if (lookahead == '-') ADVANCE(10);
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(289);
      if (lookahead == '2') ADVANCE(290);
      if (lookahead == '3') ADVANCE(291);
      if (lookahead == ':') ADVANCE(410);
      if (lookahead == '<') ADVANCE(19);
      if (lookahead == '=') ADVANCE(261);
      if (lookahead == '[') ADVANCE(20);
      if (lookahead == 'c') ADVANCE(198);
      if (lookahead == 'e') ADVANCE(244);
      if (lookahead == 'f') ADVANCE(76);
      if (lookahead == 'g') ADVANCE(115);
      if (lookahead == 'j') ADVANCE(77);
      if (lookahead == 'k') ADVANCE(213);
      if (lookahead == 'l') ADVANCE(116);
      if (lookahead == 'm') ADVANCE(117);
      if (lookahead == 'n') ADVANCE(118);
      if (lookahead == 'o') ADVANCE(143);
      if (lookahead == 'r') ADVANCE(78);
      if (lookahead == 's') ADVANCE(129);
      if (lookahead == 't') ADVANCE(205);
      if (lookahead == 'y') ADVANCE(133);
      if (lookahead == '{') ADVANCE(400);
      if (lookahead == '}') ADVANCE(401);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ') ADVANCE(3);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 4:
      if (lookahead == '"') ADVANCE(399);
      if (lookahead == '#') ADVANCE(107);
      if (lookahead == '$') ADVANCE(252);
      if (lookahead == '\'') ADVANCE(398);
      if (lookahead == '+' ||
          lookahead == '-') ADVANCE(13);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(4)
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(289);
      if (lookahead == '2') ADVANCE(290);
      if (lookahead == '3') ADVANCE(291);
      if (lookahead == 'f') ADVANCE(345);
      if (lookahead == 'g') ADVANCE(358);
      if (lookahead == 'n') ADVANCE(351);
      if (lookahead == 'o') ADVANCE(360);
      if (lookahead == 't') ADVANCE(381);
      if (lookahead == 'y') ADVANCE(355);
      if (lookahead == '}') ADVANCE(401);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 5:
      if (lookahead == '"') ADVANCE(399);
      if (lookahead == '$') ADVANCE(251);
      if (lookahead == '\'') ADVANCE(398);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(5)
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(102);
      if (lookahead == '2') ADVANCE(103);
      if (lookahead == '3') ADVANCE(104);
      if (lookahead == '=') ADVANCE(261);
      if (lookahead == '{') ADVANCE(400);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 6:
      if (lookahead == '"') ADVANCE(399);
      if (lookahead == '$') ADVANCE(252);
      if (lookahead == '\'') ADVANCE(398);
      if (lookahead == '(') ADVANCE(402);
      if (lookahead == ')') ADVANCE(403);
      if (lookahead == '+' ||
          lookahead == '-') ADVANCE(13);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(6)
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == ':') ADVANCE(410);
      if (lookahead == '[') ADVANCE(20);
      if (lookahead == 'f') ADVANCE(345);
      if (lookahead == 'n') ADVANCE(374);
      if (lookahead == 'o') ADVANCE(360);
      if (lookahead == 't') ADVANCE(381);
      if (lookahead == 'y') ADVANCE(355);
      if (lookahead == '{') ADVANCE(400);
      if (lookahead == '}') ADVANCE(401);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 7:
      if (lookahead == '"') ADVANCE(399);
      if (lookahead == '$') ADVANCE(252);
      if (lookahead == '\'') ADVANCE(398);
      if (lookahead == '+' ||
          lookahead == '-') ADVANCE(13);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(7)
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(289);
      if (lookahead == '2') ADVANCE(290);
      if (lookahead == '3') ADVANCE(291);
      if (lookahead == '=') ADVANCE(261);
      if (lookahead == 'f') ADVANCE(345);
      if (lookahead == 'n') ADVANCE(374);
      if (lookahead == 'o') ADVANCE(360);
      if (lookahead == 't') ADVANCE(381);
      if (lookahead == 'y') ADVANCE(355);
      if (lookahead == '}') ADVANCE(401);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 8:
      if (lookahead == '"') ADVANCE(399);
      if (lookahead == '\'') ADVANCE(398);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(8)
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(102);
      if (lookahead == '2') ADVANCE(103);
      if (lookahead == '3') ADVANCE(104);
      if (lookahead == 'f') ADVANCE(378);
      if (lookahead == '{') ADVANCE(400);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 9:
      if (lookahead == '#') ADVANCE(107);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(9)
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(102);
      if (lookahead == '2') ADVANCE(103);
      if (lookahead == '3') ADVANCE(104);
      if (lookahead == '[') ADVANCE(20);
      if (lookahead == 'g') ADVANCE(358);
      if (lookahead == 'n') ADVANCE(352);
      if (lookahead == '}') ADVANCE(401);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 10:
      if (lookahead == '$') ADVANCE(252);
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '>') ADVANCE(282);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ') ADVANCE(409);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 11:
      if (lookahead == '$') ADVANCE(252);
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '>') ADVANCE(282);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 12:
      if (lookahead == '$') ADVANCE(252);
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ') ADVANCE(408);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 13:
      if (lookahead == '$') ADVANCE(252);
      if (lookahead == '.') ADVANCE(249);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 14:
      if (lookahead == '$') ADVANCE(252);
      if (lookahead == '+' ||
          lookahead == '-') ADVANCE(13);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(14)
      if (lookahead == '.') ADVANCE(249);
      if (lookahead == '/') ADVANCE(15);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 15:
      if (lookahead == '*') ADVANCE(17);
      if (lookahead == '/') ADVANCE(406);
      END_STATE();
    case 16:
      if (lookahead == '*') ADVANCE(16);
      if (lookahead == '/') ADVANCE(405);
      if (lookahead != 0) ADVANCE(17);
      END_STATE();
    case 17:
      if (lookahead == '*') ADVANCE(16);
      if (lookahead != 0) ADVANCE(17);
      END_STATE();
    case 18:
      if (lookahead == '*') ADVANCE(95);
      if (lookahead == ']') ADVANCE(343);
      END_STATE();
    case 19:
      if (lookahead == '-') ADVANCE(39);
      END_STATE();
    case 20:
      if (lookahead == '/') ADVANCE(93);
      if (lookahead == 'G') ADVANCE(18);
      if (lookahead == 'J') ADVANCE(21);
      if (lookahead == 'K') ADVANCE(41);
      if (lookahead == 'P') ADVANCE(81);
      if (lookahead == 'T') ADVANCE(42);
      if (lookahead == 'V') ADVANCE(22);
      if (lookahead == 'b') ADVANCE(80);
      if (lookahead == 'c') ADVANCE(176);
      if (lookahead == 'd') ADVANCE(128);
      if (lookahead == 'e') ADVANCE(40);
      if (lookahead == 'f') ADVANCE(218);
      if (lookahead == 'g') ADVANCE(23);
      if (lookahead == 'k') ADVANCE(151);
      if (lookahead == 'm') ADVANCE(30);
      if (lookahead == 'n') ADVANCE(220);
      if (lookahead == 'p') ADVANCE(222);
      if (lookahead == 'r') ADVANCE(79);
      if (lookahead == 's') ADVANCE(43);
      if (lookahead == 'u') ADVANCE(180);
      if (lookahead == 'w') ADVANCE(119);
      END_STATE();
    case 21:
      if (lookahead == '/') ADVANCE(99);
      END_STATE();
    case 22:
      if (lookahead == '/') ADVANCE(98);
      if (lookahead == ']') ADVANCE(337);
      END_STATE();
    case 23:
      if (lookahead == '/') ADVANCE(101);
      END_STATE();
    case 24:
      if (lookahead == '/') ADVANCE(97);
      END_STATE();
    case 25:
      if (lookahead == '/') ADVANCE(155);
      END_STATE();
    case 26:
      if (lookahead == '/') ADVANCE(183);
      END_STATE();
    case 27:
      if (lookahead == '/') ADVANCE(225);
      if (lookahead == ']') ADVANCE(334);
      END_STATE();
    case 28:
      if (lookahead == '/') ADVANCE(226);
      if (lookahead == ']') ADVANCE(333);
      END_STATE();
    case 29:
      if (lookahead == '/') ADVANCE(187);
      END_STATE();
    case 30:
      if (lookahead == '2') ADVANCE(27);
      if (lookahead == ']') ADVANCE(314);
      if (lookahead == 'm') ADVANCE(47);
      if (lookahead == 'r') ADVANCE(84);
      END_STATE();
    case 31:
      if (lookahead == '2') ADVANCE(28);
      if (lookahead == ']') ADVANCE(313);
      END_STATE();
    case 32:
      if (lookahead == '2') ADVANCE(71);
      END_STATE();
    case 33:
      if (lookahead == '3') ADVANCE(52);
      END_STATE();
    case 34:
      if (lookahead == '3') ADVANCE(57);
      END_STATE();
    case 35:
      if (lookahead == '3') ADVANCE(60);
      END_STATE();
    case 36:
      if (lookahead == '3') ADVANCE(65);
      END_STATE();
    case 37:
      if (lookahead == '3') ADVANCE(67);
      END_STATE();
    case 38:
      if (lookahead == '3') ADVANCE(68);
      END_STATE();
    case 39:
      if (lookahead == '>') ADVANCE(272);
      END_STATE();
    case 40:
      if (lookahead == 'V') ADVANCE(45);
      if (lookahead == 'r') ADVANCE(152);
      END_STATE();
    case 41:
      if (lookahead == ']') ADVANCE(327);
      END_STATE();
    case 42:
      if (lookahead == ']') ADVANCE(342);
      END_STATE();
    case 43:
      if (lookahead == ']') ADVANCE(319);
      END_STATE();
    case 44:
      if (lookahead == ']') ADVANCE(328);
      END_STATE();
    case 45:
      if (lookahead == ']') ADVANCE(326);
      END_STATE();
    case 46:
      if (lookahead == ']') ADVANCE(315);
      END_STATE();
    case 47:
      if (lookahead == ']') ADVANCE(312);
      END_STATE();
    case 48:
      if (lookahead == ']') ADVANCE(317);
      END_STATE();
    case 49:
      if (lookahead == ']') ADVANCE(316);
      END_STATE();
    case 50:
      if (lookahead == ']') ADVANCE(311);
      END_STATE();
    case 51:
      if (lookahead == ']') ADVANCE(318);
      END_STATE();
    case 52:
      if (lookahead == ']') ADVANCE(320);
      END_STATE();
    case 53:
      if (lookahead == ']') ADVANCE(340);
      END_STATE();
    case 54:
      if (lookahead == ']') ADVANCE(330);
      END_STATE();
    case 55:
      if (lookahead == ']') ADVANCE(307);
      END_STATE();
    case 56:
      if (lookahead == ']') ADVANCE(308);
      END_STATE();
    case 57:
      if (lookahead == ']') ADVANCE(321);
      END_STATE();
    case 58:
      if (lookahead == ']') ADVANCE(339);
      END_STATE();
    case 59:
      if (lookahead == ']') ADVANCE(332);
      END_STATE();
    case 60:
      if (lookahead == ']') ADVANCE(324);
      END_STATE();
    case 61:
      if (lookahead == ']') ADVANCE(341);
      END_STATE();
    case 62:
      if (lookahead == ']') ADVANCE(336);
      END_STATE();
    case 63:
      if (lookahead == ']') ADVANCE(309);
      END_STATE();
    case 64:
      if (lookahead == ']') ADVANCE(310);
      END_STATE();
    case 65:
      if (lookahead == ']') ADVANCE(325);
      END_STATE();
    case 66:
      if (lookahead == ']') ADVANCE(335);
      END_STATE();
    case 67:
      if (lookahead == ']') ADVANCE(323);
      END_STATE();
    case 68:
      if (lookahead == ']') ADVANCE(322);
      END_STATE();
    case 69:
      if (lookahead == ']') ADVANCE(331);
      END_STATE();
    case 70:
      if (lookahead == ']') ADVANCE(338);
      END_STATE();
    case 71:
      if (lookahead == ']') ADVANCE(329);
      END_STATE();
    case 72:
      if (lookahead == '_') ADVANCE(142);
      END_STATE();
    case 73:
      if (lookahead == '_') ADVANCE(134);
      END_STATE();
    case 74:
      if (lookahead == '_') ADVANCE(111);
      END_STATE();
    case 75:
      if (lookahead == '_') ADVANCE(212);
      END_STATE();
    case 76:
      if (lookahead == 'a') ADVANCE(175);
      if (lookahead == 'e') ADVANCE(207);
      if (lookahead == 'o') ADVANCE(206);
      END_STATE();
    case 77:
      if (lookahead == 'a') ADVANCE(193);
      END_STATE();
    case 78:
      if (lookahead == 'a') ADVANCE(233);
      if (lookahead == 'e') ADVANCE(88);
      END_STATE();
    case 79:
      if (lookahead == 'a') ADVANCE(109);
      END_STATE();
    case 80:
      if (lookahead == 'a') ADVANCE(211);
      END_STATE();
    case 81:
      if (lookahead == 'a') ADVANCE(44);
      END_STATE();
    case 82:
      if (lookahead == 'a') ADVANCE(174);
      END_STATE();
    case 83:
      if (lookahead == 'a') ADVANCE(234);
      END_STATE();
    case 84:
      if (lookahead == 'a') ADVANCE(113);
      END_STATE();
    case 85:
      if (lookahead == 'a') ADVANCE(235);
      END_STATE();
    case 86:
      if (lookahead == 'a') ADVANCE(114);
      END_STATE();
    case 87:
      if (lookahead == 'a') ADVANCE(237);
      END_STATE();
    case 88:
      if (lookahead == 'a') ADVANCE(100);
      if (lookahead == 'f') ADVANCE(74);
      END_STATE();
    case 89:
      if (lookahead == 'b') ADVANCE(273);
      END_STATE();
    case 90:
      if (lookahead == 'b') ADVANCE(139);
      END_STATE();
    case 91:
      if (lookahead == 'c') ADVANCE(162);
      END_STATE();
    case 92:
      if (lookahead == 'c') ADVANCE(276);
      END_STATE();
    case 93:
      if (lookahead == 'c') ADVANCE(184);
      if (lookahead == 'm') ADVANCE(33);
      END_STATE();
    case 94:
      if (lookahead == 'c') ADVANCE(232);
      END_STATE();
    case 95:
      if (lookahead == 'c') ADVANCE(181);
      END_STATE();
    case 96:
      if (lookahead == 'c') ADVANCE(170);
      END_STATE();
    case 97:
      if (lookahead == 'c') ADVANCE(177);
      END_STATE();
    case 98:
      if (lookahead == 'c') ADVANCE(182);
      if (lookahead == 'm') ADVANCE(53);
      END_STATE();
    case 99:
      if (lookahead == 'c') ADVANCE(185);
      if (lookahead == 'k') ADVANCE(154);
      if (lookahead == 'm') ADVANCE(35);
      END_STATE();
    case 100:
      if (lookahead == 'c') ADVANCE(236);
      END_STATE();
    case 101:
      if (lookahead == 'c') ADVANCE(186);
      END_STATE();
    case 102:
      if (lookahead == 'd') ADVANCE(395);
      END_STATE();
    case 103:
      if (lookahead == 'd') ADVANCE(396);
      END_STATE();
    case 104:
      if (lookahead == 'd') ADVANCE(397);
      END_STATE();
    case 105:
      if (lookahead == 'd') ADVANCE(158);
      END_STATE();
    case 106:
      if (lookahead == 'd') ADVANCE(132);
      if (lookahead == 'e') ADVANCE(171);
      if (lookahead == 'i') ADVANCE(144);
      END_STATE();
    case 107:
      if (lookahead == 'd') ADVANCE(132);
      if (lookahead == 'i') ADVANCE(144);
      END_STATE();
    case 108:
      if (lookahead == 'd') ADVANCE(135);
      if (lookahead == 'n') ADVANCE(110);
      END_STATE();
    case 109:
      if (lookahead == 'd') ADVANCE(56);
      END_STATE();
    case 110:
      if (lookahead == 'd') ADVANCE(138);
      END_STATE();
    case 111:
      if (lookahead == 'd') ADVANCE(130);
      END_STATE();
    case 112:
      if (lookahead == 'd') ADVANCE(125);
      END_STATE();
    case 113:
      if (lookahead == 'd') ADVANCE(63);
      END_STATE();
    case 114:
      if (lookahead == 'd') ADVANCE(64);
      END_STATE();
    case 115:
      if (lookahead == 'e') ADVANCE(192);
      END_STATE();
    case 116:
      if (lookahead == 'e') ADVANCE(242);
      END_STATE();
    case 117:
      if (lookahead == 'e') ADVANCE(230);
      END_STATE();
    case 118:
      if (lookahead == 'e') ADVANCE(243);
      if (lookahead == 'o') ADVANCE(301);
      END_STATE();
    case 119:
      if (lookahead == 'e') ADVANCE(90);
      END_STATE();
    case 120:
      if (lookahead == 'e') ADVANCE(411);
      END_STATE();
    case 121:
      if (lookahead == 'e') ADVANCE(295);
      END_STATE();
    case 122:
      if (lookahead == 'e') ADVANCE(286);
      END_STATE();
    case 123:
      if (lookahead == 'e') ADVANCE(297);
      END_STATE();
    case 124:
      if (lookahead == 'e') ADVANCE(257);
      END_STATE();
    case 125:
      if (lookahead == 'e') ADVANCE(256);
      END_STATE();
    case 126:
      if (lookahead == 'e') ADVANCE(266);
      END_STATE();
    case 127:
      if (lookahead == 'e') ADVANCE(412);
      END_STATE();
    case 128:
      if (lookahead == 'e') ADVANCE(153);
      if (lookahead == 'y') ADVANCE(195);
      END_STATE();
    case 129:
      if (lookahead == 'e') ADVANCE(94);
      END_STATE();
    case 130:
      if (lookahead == 'e') ADVANCE(194);
      END_STATE();
    case 131:
      if (lookahead == 'e') ADVANCE(241);
      END_STATE();
    case 132:
      if (lookahead == 'e') ADVANCE(149);
      END_STATE();
    case 133:
      if (lookahead == 'e') ADVANCE(214);
      END_STATE();
    case 134:
      if (lookahead == 'e') ADVANCE(240);
      END_STATE();
    case 135:
      if (lookahead == 'e') ADVANCE(147);
      END_STATE();
    case 136:
      if (lookahead == 'e') ADVANCE(209);
      END_STATE();
    case 137:
      if (lookahead == 'e') ADVANCE(167);
      END_STATE();
    case 138:
      if (lookahead == 'e') ADVANCE(148);
      END_STATE();
    case 139:
      if (lookahead == 'e') ADVANCE(210);
      END_STATE();
    case 140:
      if (lookahead == 'e') ADVANCE(208);
      END_STATE();
    case 141:
      if (lookahead == 'e') ADVANCE(227);
      END_STATE();
    case 142:
      if (lookahead == 'e') ADVANCE(197);
      END_STATE();
    case 143:
      if (lookahead == 'f') ADVANCE(145);
      if (lookahead == 'n') ADVANCE(303);
      END_STATE();
    case 144:
      if (lookahead == 'f') ADVANCE(108);
      if (lookahead == 'n') ADVANCE(96);
      END_STATE();
    case 145:
      if (lookahead == 'f') ADVANCE(305);
      END_STATE();
    case 146:
      if (lookahead == 'f') ADVANCE(260);
      END_STATE();
    case 147:
      if (lookahead == 'f') ADVANCE(258);
      END_STATE();
    case 148:
      if (lookahead == 'f') ADVANCE(259);
      END_STATE();
    case 149:
      if (lookahead == 'f') ADVANCE(164);
      END_STATE();
    case 150:
      if (lookahead == 'g') ADVANCE(246);
      END_STATE();
    case 151:
      if (lookahead == 'g') ADVANCE(29);
      END_STATE();
    case 152:
      if (lookahead == 'g') ADVANCE(216);
      END_STATE();
    case 153:
      if (lookahead == 'g') ADVANCE(55);
      END_STATE();
    case 154:
      if (lookahead == 'g') ADVANCE(59);
      END_STATE();
    case 155:
      if (lookahead == 'g') ADVANCE(69);
      END_STATE();
    case 156:
      if (lookahead == 'i') ADVANCE(72);
      END_STATE();
    case 157:
      if (lookahead == 'i') ADVANCE(199);
      END_STATE();
    case 158:
      if (lookahead == 'i') ADVANCE(146);
      END_STATE();
    case 159:
      if (lookahead == 'i') ADVANCE(229);
      END_STATE();
    case 160:
      if (lookahead == 'i') ADVANCE(92);
      END_STATE();
    case 161:
      if (lookahead == 'i') ADVANCE(228);
      END_STATE();
    case 162:
      if (lookahead == 'i') ADVANCE(231);
      END_STATE();
    case 163:
      if (lookahead == 'i') ADVANCE(200);
      END_STATE();
    case 164:
      if (lookahead == 'i') ADVANCE(196);
      END_STATE();
    case 165:
      if (lookahead == 'i') ADVANCE(203);
      END_STATE();
    case 166:
      if (lookahead == 'i') ADVANCE(204);
      END_STATE();
    case 167:
      if (lookahead == 'l') ADVANCE(283);
      END_STATE();
    case 168:
      if (lookahead == 'l') ADVANCE(202);
      END_STATE();
    case 169:
      if (lookahead == 'l') ADVANCE(172);
      if (lookahead == 'u') ADVANCE(168);
      END_STATE();
    case 170:
      if (lookahead == 'l') ADVANCE(239);
      END_STATE();
    case 171:
      if (lookahead == 'l') ADVANCE(223);
      if (lookahead == 'n') ADVANCE(105);
      END_STATE();
    case 172:
      if (lookahead == 'l') ADVANCE(161);
      END_STATE();
    case 173:
      if (lookahead == 'l') ADVANCE(160);
      END_STATE();
    case 174:
      if (lookahead == 'l') ADVANCE(173);
      END_STATE();
    case 175:
      if (lookahead == 'l') ADVANCE(224);
      END_STATE();
    case 176:
      if (lookahead == 'm') ADVANCE(31);
      END_STATE();
    case 177:
      if (lookahead == 'm') ADVANCE(32);
      END_STATE();
    case 178:
      if (lookahead == 'm') ADVANCE(89);
      END_STATE();
    case 179:
      if (lookahead == 'm') ADVANCE(156);
      END_STATE();
    case 180:
      if (lookahead == 'm') ADVANCE(50);
      if (lookahead == 'r') ADVANCE(86);
      if (lookahead == 's') ADVANCE(51);
      END_STATE();
    case 181:
      if (lookahead == 'm') ADVANCE(58);
      END_STATE();
    case 182:
      if (lookahead == 'm') ADVANCE(61);
      END_STATE();
    case 183:
      if (lookahead == 'm') ADVANCE(70);
      END_STATE();
    case 184:
      if (lookahead == 'm') ADVANCE(34);
      END_STATE();
    case 185:
      if (lookahead == 'm') ADVANCE(36);
      END_STATE();
    case 186:
      if (lookahead == 'm') ADVANCE(37);
      END_STATE();
    case 187:
      if (lookahead == 'm') ADVANCE(38);
      END_STATE();
    case 188:
      if (lookahead == 'n') ADVANCE(275);
      END_STATE();
    case 189:
      if (lookahead == 'n') ADVANCE(268);
      END_STATE();
    case 190:
      if (lookahead == 'n') ADVANCE(270);
      END_STATE();
    case 191:
      if (lookahead == 'n') ADVANCE(280);
      END_STATE();
    case 192:
      if (lookahead == 'n') ADVANCE(136);
      if (lookahead == 't') ADVANCE(284);
      END_STATE();
    case 193:
      if (lookahead == 'n') ADVANCE(131);
      END_STATE();
    case 194:
      if (lookahead == 'n') ADVANCE(221);
      END_STATE();
    case 195:
      if (lookahead == 'n') ADVANCE(141);
      END_STATE();
    case 196:
      if (lookahead == 'n') ADVANCE(124);
      END_STATE();
    case 197:
      if (lookahead == 'n') ADVANCE(140);
      END_STATE();
    case 198:
      if (lookahead == 'o') ADVANCE(169);
      if (lookahead == 'r') ADVANCE(201);
      END_STATE();
    case 199:
      if (lookahead == 'o') ADVANCE(188);
      END_STATE();
    case 200:
      if (lookahead == 'o') ADVANCE(189);
      END_STATE();
    case 201:
      if (lookahead == 'o') ADVANCE(219);
      END_STATE();
    case 202:
      if (lookahead == 'o') ADVANCE(178);
      END_STATE();
    case 203:
      if (lookahead == 'o') ADVANCE(190);
      END_STATE();
    case 204:
      if (lookahead == 'o') ADVANCE(191);
      END_STATE();
    case 205:
      if (lookahead == 'r') ADVANCE(238);
      END_STATE();
    case 206:
      if (lookahead == 'r') ADVANCE(264);
      END_STATE();
    case 207:
      if (lookahead == 'r') ADVANCE(179);
      END_STATE();
    case 208:
      if (lookahead == 'r') ADVANCE(150);
      END_STATE();
    case 209:
      if (lookahead == 'r') ADVANCE(83);
      END_STATE();
    case 210:
      if (lookahead == 'r') ADVANCE(217);
      END_STATE();
    case 211:
      if (lookahead == 'r') ADVANCE(54);
      END_STATE();
    case 212:
      if (lookahead == 'r') ADVANCE(85);
      END_STATE();
    case 213:
      if (lookahead == 's') ADVANCE(277);
      END_STATE();
    case 214:
      if (lookahead == 's') ADVANCE(299);
      END_STATE();
    case 215:
      if (lookahead == 's') ADVANCE(274);
      END_STATE();
    case 216:
      if (lookahead == 's') ADVANCE(25);
      END_STATE();
    case 217:
      if (lookahead == 's') ADVANCE(26);
      END_STATE();
    case 218:
      if (lookahead == 's') ADVANCE(46);
      END_STATE();
    case 219:
      if (lookahead == 's') ADVANCE(215);
      END_STATE();
    case 220:
      if (lookahead == 's') ADVANCE(48);
      END_STATE();
    case 221:
      if (lookahead == 's') ADVANCE(159);
      END_STATE();
    case 222:
      if (lookahead == 's') ADVANCE(49);
      END_STATE();
    case 223:
      if (lookahead == 's') ADVANCE(122);
      END_STATE();
    case 224:
      if (lookahead == 's') ADVANCE(123);
      END_STATE();
    case 225:
      if (lookahead == 's') ADVANCE(62);
      END_STATE();
    case 226:
      if (lookahead == 's') ADVANCE(66);
      END_STATE();
    case 227:
      if (lookahead == 's') ADVANCE(24);
      END_STATE();
    case 228:
      if (lookahead == 's') ADVANCE(165);
      END_STATE();
    case 229:
      if (lookahead == 't') ADVANCE(245);
      END_STATE();
    case 230:
      if (lookahead == 't') ADVANCE(82);
      END_STATE();
    case 231:
      if (lookahead == 't') ADVANCE(87);
      END_STATE();
    case 232:
      if (lookahead == 't') ADVANCE(157);
      END_STATE();
    case 233:
      if (lookahead == 't') ADVANCE(120);
      END_STATE();
    case 234:
      if (lookahead == 't') ADVANCE(126);
      END_STATE();
    case 235:
      if (lookahead == 't') ADVANCE(127);
      END_STATE();
    case 236:
      if (lookahead == 't') ADVANCE(163);
      END_STATE();
    case 237:
      if (lookahead == 't') ADVANCE(166);
      END_STATE();
    case 238:
      if (lookahead == 'u') ADVANCE(121);
      END_STATE();
    case 239:
      if (lookahead == 'u') ADVANCE(112);
      END_STATE();
    case 240:
      if (lookahead == 'v') ADVANCE(278);
      END_STATE();
    case 241:
      if (lookahead == 'v') ADVANCE(75);
      END_STATE();
    case 242:
      if (lookahead == 'v') ADVANCE(137);
      END_STATE();
    case 243:
      if (lookahead == 'w') ADVANCE(262);
      END_STATE();
    case 244:
      if (lookahead == 'x') ADVANCE(91);
      END_STATE();
    case 245:
      if (lookahead == 'y') ADVANCE(279);
      END_STATE();
    case 246:
      if (lookahead == 'y') ADVANCE(73);
      END_STATE();
    case 247:
      if (lookahead == '+' ||
          lookahead == '-') ADVANCE(250);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(294);
      END_STATE();
    case 248:
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(248)
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(102);
      if (lookahead == '2') ADVANCE(103);
      if (lookahead == '3') ADVANCE(104);
      if (lookahead == 'c') ADVANCE(376);
      if (lookahead == 'e') ADVANCE(393);
      if (lookahead == 'r') ADVANCE(359);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 249:
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(293);
      END_STATE();
    case 250:
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(294);
      END_STATE();
    case 251:
      if (!sym_define_key_character_set_1(lookahead)) ADVANCE(287);
      END_STATE();
    case 252:
      if (!sym_define_key_character_set_1(lookahead)) ADVANCE(288);
      END_STATE();
    case 253:
      if (lookahead != 0 &&
          lookahead != '\r') ADVANCE(406);
      if (lookahead == '\r') ADVANCE(407);
      END_STATE();
    case 254:
      if (eof) ADVANCE(255);
      if (lookahead == '#') ADVANCE(106);
      if (lookahead == '$') ADVANCE(251);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ' ||
          lookahead == ',') SKIP(254)
      if (lookahead == '/') ADVANCE(15);
      if (lookahead == '1') ADVANCE(102);
      if (lookahead == '2') ADVANCE(103);
      if (lookahead == '3') ADVANCE(104);
      if (lookahead == '[') ADVANCE(20);
      if (lookahead == 'g') ADVANCE(350);
      if (lookahead == 'n') ADVANCE(352);
      if (lookahead == '{') ADVANCE(400);
      if (('A' <= lookahead && lookahead <= 'Z') ||
          lookahead == '_' ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 255:
      ACCEPT_TOKEN(ts_builtin_sym_end);
      END_STATE();
    case 256:
      ACCEPT_TOKEN(anon_sym_POUNDinclude);
      END_STATE();
    case 257:
      ACCEPT_TOKEN(anon_sym_POUNDdefine);
      END_STATE();
    case 258:
      ACCEPT_TOKEN(anon_sym_POUNDifdef);
      END_STATE();
    case 259:
      ACCEPT_TOKEN(anon_sym_POUNDifndef);
      END_STATE();
    case 260:
      ACCEPT_TOKEN(anon_sym_POUNDendif);
      END_STATE();
    case 261:
      ACCEPT_TOKEN(anon_sym_EQ);
      END_STATE();
    case 262:
      ACCEPT_TOKEN(anon_sym_new);
      END_STATE();
    case 263:
      ACCEPT_TOKEN(anon_sym_new);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 264:
      ACCEPT_TOKEN(anon_sym_for);
      END_STATE();
    case 265:
      ACCEPT_TOKEN(anon_sym_for);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 266:
      ACCEPT_TOKEN(anon_sym_generate);
      END_STATE();
    case 267:
      ACCEPT_TOKEN(anon_sym_generate);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 268:
      ACCEPT_TOKEN(anon_sym_reaction);
      END_STATE();
    case 269:
      ACCEPT_TOKEN(anon_sym_reaction);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 270:
      ACCEPT_TOKEN(anon_sym_collision);
      END_STATE();
    case 271:
      ACCEPT_TOKEN(anon_sym_collision);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 272:
      ACCEPT_TOKEN(anon_sym_LT_DASH_GT);
      END_STATE();
    case 273:
      ACCEPT_TOKEN(anon_sym_coulomb);
      END_STATE();
    case 274:
      ACCEPT_TOKEN(anon_sym_cross);
      END_STATE();
    case 275:
      ACCEPT_TOKEN(anon_sym_section);
      END_STATE();
    case 276:
      ACCEPT_TOKEN(anon_sym_metallic);
      END_STATE();
    case 277:
      ACCEPT_TOKEN(anon_sym_ks);
      END_STATE();
    case 278:
      ACCEPT_TOKEN(anon_sym_fermi_energy_ev);
      END_STATE();
    case 279:
      ACCEPT_TOKEN(anon_sym_ref_density);
      END_STATE();
    case 280:
      ACCEPT_TOKEN(anon_sym_excitation);
      END_STATE();
    case 281:
      ACCEPT_TOKEN(anon_sym_excitation);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 282:
      ACCEPT_TOKEN(anon_sym_DASH_GT);
      END_STATE();
    case 283:
      ACCEPT_TOKEN(anon_sym_level);
      END_STATE();
    case 284:
      ACCEPT_TOKEN(anon_sym_get);
      END_STATE();
    case 285:
      ACCEPT_TOKEN(anon_sym_get);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 286:
      ACCEPT_TOKEN(anon_sym_POUNDelse);
      END_STATE();
    case 287:
      ACCEPT_TOKEN(sym_define_key);
      if (!sym_define_key_character_set_1(lookahead)) ADVANCE(287);
      END_STATE();
    case 288:
      ACCEPT_TOKEN(sym_define_ref);
      if (!sym_define_key_character_set_1(lookahead)) ADVANCE(288);
      END_STATE();
    case 289:
      ACCEPT_TOKEN(sym_decimal);
      if (lookahead == '.') ADVANCE(293);
      if (lookahead == 'd') ADVANCE(395);
      if (lookahead == 'E' ||
          lookahead == 'e') ADVANCE(247);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 290:
      ACCEPT_TOKEN(sym_decimal);
      if (lookahead == '.') ADVANCE(293);
      if (lookahead == 'd') ADVANCE(396);
      if (lookahead == 'E' ||
          lookahead == 'e') ADVANCE(247);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 291:
      ACCEPT_TOKEN(sym_decimal);
      if (lookahead == '.') ADVANCE(293);
      if (lookahead == 'd') ADVANCE(397);
      if (lookahead == 'E' ||
          lookahead == 'e') ADVANCE(247);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 292:
      ACCEPT_TOKEN(sym_decimal);
      if (lookahead == '.') ADVANCE(293);
      if (lookahead == 'E' ||
          lookahead == 'e') ADVANCE(247);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(292);
      END_STATE();
    case 293:
      ACCEPT_TOKEN(sym_decimal);
      if (lookahead == 'E' ||
          lookahead == 'e') ADVANCE(247);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(293);
      END_STATE();
    case 294:
      ACCEPT_TOKEN(sym_decimal);
      if (('0' <= lookahead && lookahead <= '9')) ADVANCE(294);
      END_STATE();
    case 295:
      ACCEPT_TOKEN(anon_sym_true);
      END_STATE();
    case 296:
      ACCEPT_TOKEN(anon_sym_true);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 297:
      ACCEPT_TOKEN(anon_sym_false);
      END_STATE();
    case 298:
      ACCEPT_TOKEN(anon_sym_false);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 299:
      ACCEPT_TOKEN(anon_sym_yes);
      END_STATE();
    case 300:
      ACCEPT_TOKEN(anon_sym_yes);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 301:
      ACCEPT_TOKEN(anon_sym_no);
      END_STATE();
    case 302:
      ACCEPT_TOKEN(anon_sym_no);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 303:
      ACCEPT_TOKEN(anon_sym_on);
      END_STATE();
    case 304:
      ACCEPT_TOKEN(anon_sym_on);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 305:
      ACCEPT_TOKEN(anon_sym_off);
      END_STATE();
    case 306:
      ACCEPT_TOKEN(anon_sym_off);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 307:
      ACCEPT_TOKEN(anon_sym_LBRACKdeg_RBRACK);
      END_STATE();
    case 308:
      ACCEPT_TOKEN(anon_sym_LBRACKrad_RBRACK);
      END_STATE();
    case 309:
      ACCEPT_TOKEN(anon_sym_LBRACKmrad_RBRACK);
      END_STATE();
    case 310:
      ACCEPT_TOKEN(anon_sym_LBRACKurad_RBRACK);
      END_STATE();
    case 311:
      ACCEPT_TOKEN(anon_sym_LBRACKum_RBRACK);
      END_STATE();
    case 312:
      ACCEPT_TOKEN(anon_sym_LBRACKmm_RBRACK);
      END_STATE();
    case 313:
      ACCEPT_TOKEN(anon_sym_LBRACKcm_RBRACK);
      END_STATE();
    case 314:
      ACCEPT_TOKEN(anon_sym_LBRACKm_RBRACK);
      END_STATE();
    case 315:
      ACCEPT_TOKEN(anon_sym_LBRACKfs_RBRACK);
      END_STATE();
    case 316:
      ACCEPT_TOKEN(anon_sym_LBRACKps_RBRACK);
      END_STATE();
    case 317:
      ACCEPT_TOKEN(anon_sym_LBRACKns_RBRACK);
      END_STATE();
    case 318:
      ACCEPT_TOKEN(anon_sym_LBRACKus_RBRACK);
      END_STATE();
    case 319:
      ACCEPT_TOKEN(anon_sym_LBRACKs_RBRACK);
      END_STATE();
    case 320:
      ACCEPT_TOKEN(anon_sym_LBRACK_SLASHm3_RBRACK);
      END_STATE();
    case 321:
      ACCEPT_TOKEN(anon_sym_LBRACK_SLASHcm3_RBRACK);
      END_STATE();
    case 322:
      ACCEPT_TOKEN(anon_sym_LBRACKkg_SLASHm3_RBRACK);
      END_STATE();
    case 323:
      ACCEPT_TOKEN(anon_sym_LBRACKg_SLASHcm3_RBRACK);
      END_STATE();
    case 324:
      ACCEPT_TOKEN(anon_sym_LBRACKJ_SLASHm3_RBRACK);
      END_STATE();
    case 325:
      ACCEPT_TOKEN(anon_sym_LBRACKJ_SLASHcm3_RBRACK);
      END_STATE();
    case 326:
      ACCEPT_TOKEN(anon_sym_LBRACKeV_RBRACK);
      END_STATE();
    case 327:
      ACCEPT_TOKEN(anon_sym_LBRACKK_RBRACK);
      END_STATE();
    case 328:
      ACCEPT_TOKEN(anon_sym_LBRACKPa_RBRACK);
      END_STATE();
    case 329:
      ACCEPT_TOKEN(anon_sym_LBRACKdynes_SLASHcm2_RBRACK);
      END_STATE();
    case 330:
      ACCEPT_TOKEN(anon_sym_LBRACKbar_RBRACK);
      END_STATE();
    case 331:
      ACCEPT_TOKEN(anon_sym_LBRACKergs_SLASHg_RBRACK);
      END_STATE();
    case 332:
      ACCEPT_TOKEN(anon_sym_LBRACKJ_SLASHkg_RBRACK);
      END_STATE();
    case 333:
      ACCEPT_TOKEN(anon_sym_LBRACKcm2_RBRACK);
      END_STATE();
    case 334:
      ACCEPT_TOKEN(anon_sym_LBRACKm2_RBRACK);
      END_STATE();
    case 335:
      ACCEPT_TOKEN(anon_sym_LBRACKcm2_SLASHs_RBRACK);
      END_STATE();
    case 336:
      ACCEPT_TOKEN(anon_sym_LBRACKm2_SLASHs_RBRACK);
      END_STATE();
    case 337:
      ACCEPT_TOKEN(anon_sym_LBRACKV_RBRACK);
      END_STATE();
    case 338:
      ACCEPT_TOKEN(anon_sym_LBRACKwebers_SLASHm_RBRACK);
      END_STATE();
    case 339:
      ACCEPT_TOKEN(anon_sym_LBRACKG_STARcm_RBRACK);
      END_STATE();
    case 340:
      ACCEPT_TOKEN(anon_sym_LBRACKV_SLASHm_RBRACK);
      END_STATE();
    case 341:
      ACCEPT_TOKEN(anon_sym_LBRACKV_SLASHcm_RBRACK);
      END_STATE();
    case 342:
      ACCEPT_TOKEN(anon_sym_LBRACKT_RBRACK);
      END_STATE();
    case 343:
      ACCEPT_TOKEN(anon_sym_LBRACKG_RBRACK);
      END_STATE();
    case 344:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'a') ADVANCE(388);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('b' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 345:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'a') ADVANCE(367);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('b' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 346:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'a') ADVANCE(349);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('b' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 347:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'a') ADVANCE(390);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('b' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 348:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'c') ADVANCE(362);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 349:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'c') ADVANCE(387);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 350:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(373);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 351:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(392);
      if (lookahead == 'o') ADVANCE(302);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 352:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(392);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 353:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(380);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 354:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(267);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 355:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(383);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 356:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(296);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 357:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(298);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 358:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(386);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 359:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'e') ADVANCE(346);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 360:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'f') ADVANCE(361);
      if (lookahead == 'n') ADVANCE(304);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 361:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'f') ADVANCE(306);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 362:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'i') ADVANCE(389);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 363:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'i') ADVANCE(385);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 364:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'i') ADVANCE(375);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 365:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'i') ADVANCE(377);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 366:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'i') ADVANCE(379);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 367:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'l') ADVANCE(384);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 368:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'l') ADVANCE(363);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 369:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'l') ADVANCE(368);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 370:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'n') ADVANCE(269);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 371:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'n') ADVANCE(271);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 372:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'n') ADVANCE(281);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 373:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'n') ADVANCE(353);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 374:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'o') ADVANCE(302);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 375:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'o') ADVANCE(370);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 376:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'o') ADVANCE(369);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 377:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'o') ADVANCE(371);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 378:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'o') ADVANCE(382);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 379:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'o') ADVANCE(372);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 380:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'r') ADVANCE(344);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 381:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'r') ADVANCE(391);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 382:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'r') ADVANCE(265);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 383:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 's') ADVANCE(300);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 384:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 's') ADVANCE(357);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 385:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 's') ADVANCE(365);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 386:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 't') ADVANCE(285);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 387:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 't') ADVANCE(364);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 388:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 't') ADVANCE(354);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 389:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 't') ADVANCE(347);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 390:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 't') ADVANCE(366);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 391:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'u') ADVANCE(356);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 392:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'w') ADVANCE(263);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 393:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == 'x') ADVANCE(348);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 394:
      ACCEPT_TOKEN(sym_identifier);
      if (lookahead == '+' ||
          lookahead == '-' ||
          lookahead == '.' ||
          ('0' <= lookahead && lookahead <= '9') ||
          ('A' <= lookahead && lookahead <= '[') ||
          (']' <= lookahead && lookahead <= '_') ||
          ('a' <= lookahead && lookahead <= 'z')) ADVANCE(394);
      END_STATE();
    case 395:
      ACCEPT_TOKEN(anon_sym_1d);
      END_STATE();
    case 396:
      ACCEPT_TOKEN(anon_sym_2d);
      END_STATE();
    case 397:
      ACCEPT_TOKEN(anon_sym_3d);
      END_STATE();
    case 398:
      ACCEPT_TOKEN(anon_sym_SQUOTE);
      END_STATE();
    case 399:
      ACCEPT_TOKEN(anon_sym_DQUOTE);
      END_STATE();
    case 400:
      ACCEPT_TOKEN(anon_sym_LBRACE);
      END_STATE();
    case 401:
      ACCEPT_TOKEN(anon_sym_RBRACE);
      END_STATE();
    case 402:
      ACCEPT_TOKEN(anon_sym_LPAREN);
      END_STATE();
    case 403:
      ACCEPT_TOKEN(anon_sym_RPAREN);
      END_STATE();
    case 404:
      ACCEPT_TOKEN(aux_sym__define_value_token1);
      if (lookahead == '\n') ADVANCE(404);
      if (lookahead == '\r') ADVANCE(1);
      END_STATE();
    case 405:
      ACCEPT_TOKEN(sym_comment);
      END_STATE();
    case 406:
      ACCEPT_TOKEN(sym_comment);
      if (lookahead == '\\') ADVANCE(253);
      if (lookahead != 0 &&
          lookahead != '\n') ADVANCE(406);
      END_STATE();
    case 407:
      ACCEPT_TOKEN(sym_comment);
      if (lookahead != 0 &&
          lookahead != '\\') ADVANCE(406);
      if (lookahead == '\\') ADVANCE(253);
      END_STATE();
    case 408:
      ACCEPT_TOKEN(aux_sym__pterm_token1);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ') ADVANCE(408);
      END_STATE();
    case 409:
      ACCEPT_TOKEN(aux_sym__nterm_token1);
      if (lookahead == '\t' ||
          lookahead == '\n' ||
          lookahead == '\r' ||
          lookahead == ' ') ADVANCE(409);
      END_STATE();
    case 410:
      ACCEPT_TOKEN(anon_sym_COLON);
      END_STATE();
    case 411:
      ACCEPT_TOKEN(anon_sym_rate);
      END_STATE();
    case 412:
      ACCEPT_TOKEN(anon_sym_janev_rate);
      END_STATE();
    default:
      return false;
  }
}

static const TSLexMode ts_lex_modes[STATE_COUNT] = {
  [0] = {.lex_state = 0},
  [1] = {.lex_state = 254},
  [2] = {.lex_state = 1},
  [3] = {.lex_state = 6},
  [4] = {.lex_state = 254},
  [5] = {.lex_state = 9},
  [6] = {.lex_state = 1},
  [7] = {.lex_state = 1},
  [8] = {.lex_state = 1},
  [9] = {.lex_state = 4},
  [10] = {.lex_state = 1},
  [11] = {.lex_state = 1},
  [12] = {.lex_state = 254},
  [13] = {.lex_state = 254},
  [14] = {.lex_state = 254},
  [15] = {.lex_state = 254},
  [16] = {.lex_state = 254},
  [17] = {.lex_state = 254},
  [18] = {.lex_state = 254},
  [19] = {.lex_state = 254},
  [20] = {.lex_state = 9},
  [21] = {.lex_state = 9},
  [22] = {.lex_state = 9},
  [23] = {.lex_state = 9},
  [24] = {.lex_state = 9},
  [25] = {.lex_state = 9},
  [26] = {.lex_state = 9},
  [27] = {.lex_state = 1},
  [28] = {.lex_state = 1},
  [29] = {.lex_state = 6},
  [30] = {.lex_state = 6},
  [31] = {.lex_state = 6},
  [32] = {.lex_state = 1},
  [33] = {.lex_state = 1},
  [34] = {.lex_state = 1},
  [35] = {.lex_state = 6},
  [36] = {.lex_state = 1},
  [37] = {.lex_state = 1},
  [38] = {.lex_state = 1},
  [39] = {.lex_state = 1},
  [40] = {.lex_state = 1},
  [41] = {.lex_state = 1},
  [42] = {.lex_state = 1},
  [43] = {.lex_state = 1},
  [44] = {.lex_state = 6},
  [45] = {.lex_state = 6},
  [46] = {.lex_state = 1},
  [47] = {.lex_state = 6},
  [48] = {.lex_state = 1},
  [49] = {.lex_state = 1},
  [50] = {.lex_state = 1},
  [51] = {.lex_state = 1},
  [52] = {.lex_state = 1},
  [53] = {.lex_state = 1},
  [54] = {.lex_state = 1},
  [55] = {.lex_state = 1},
  [56] = {.lex_state = 1},
  [57] = {.lex_state = 1},
  [58] = {.lex_state = 1},
  [59] = {.lex_state = 6},
  [60] = {.lex_state = 6},
  [61] = {.lex_state = 1},
  [62] = {.lex_state = 1},
  [63] = {.lex_state = 1},
  [64] = {.lex_state = 1},
  [65] = {.lex_state = 6},
  [66] = {.lex_state = 1},
  [67] = {.lex_state = 1},
  [68] = {.lex_state = 1},
  [69] = {.lex_state = 1},
  [70] = {.lex_state = 6},
  [71] = {.lex_state = 6},
  [72] = {.lex_state = 6},
  [73] = {.lex_state = 6},
  [74] = {.lex_state = 6},
  [75] = {.lex_state = 6},
  [76] = {.lex_state = 7},
  [77] = {.lex_state = 254},
  [78] = {.lex_state = 254},
  [79] = {.lex_state = 254},
  [80] = {.lex_state = 254},
  [81] = {.lex_state = 254},
  [82] = {.lex_state = 254},
  [83] = {.lex_state = 254},
  [84] = {.lex_state = 254},
  [85] = {.lex_state = 254},
  [86] = {.lex_state = 254},
  [87] = {.lex_state = 254},
  [88] = {.lex_state = 254},
  [89] = {.lex_state = 254},
  [90] = {.lex_state = 254},
  [91] = {.lex_state = 254},
  [92] = {.lex_state = 254},
  [93] = {.lex_state = 254},
  [94] = {.lex_state = 6},
  [95] = {.lex_state = 6},
  [96] = {.lex_state = 254},
  [97] = {.lex_state = 254},
  [98] = {.lex_state = 6},
  [99] = {.lex_state = 254},
  [100] = {.lex_state = 254},
  [101] = {.lex_state = 254},
  [102] = {.lex_state = 254},
  [103] = {.lex_state = 254},
  [104] = {.lex_state = 254},
  [105] = {.lex_state = 254},
  [106] = {.lex_state = 254},
  [107] = {.lex_state = 254},
  [108] = {.lex_state = 254},
  [109] = {.lex_state = 6},
  [110] = {.lex_state = 6},
  [111] = {.lex_state = 254},
  [112] = {.lex_state = 254},
  [113] = {.lex_state = 254},
  [114] = {.lex_state = 6},
  [115] = {.lex_state = 254},
  [116] = {.lex_state = 254},
  [117] = {.lex_state = 9},
  [118] = {.lex_state = 9},
  [119] = {.lex_state = 9},
  [120] = {.lex_state = 9},
  [121] = {.lex_state = 9},
  [122] = {.lex_state = 9},
  [123] = {.lex_state = 9},
  [124] = {.lex_state = 9},
  [125] = {.lex_state = 9},
  [126] = {.lex_state = 9},
  [127] = {.lex_state = 9},
  [128] = {.lex_state = 9},
  [129] = {.lex_state = 9},
  [130] = {.lex_state = 9},
  [131] = {.lex_state = 9},
  [132] = {.lex_state = 9},
  [133] = {.lex_state = 9},
  [134] = {.lex_state = 9},
  [135] = {.lex_state = 9},
  [136] = {.lex_state = 9},
  [137] = {.lex_state = 8},
  [138] = {.lex_state = 248},
  [139] = {.lex_state = 5},
  [140] = {.lex_state = 8},
  [141] = {.lex_state = 5},
  [142] = {.lex_state = 248},
  [143] = {.lex_state = 0},
  [144] = {.lex_state = 5},
  [145] = {.lex_state = 0},
  [146] = {.lex_state = 0},
  [147] = {.lex_state = 0},
  [148] = {.lex_state = 0},
  [149] = {.lex_state = 8},
  [150] = {.lex_state = 5},
  [151] = {.lex_state = 5},
  [152] = {.lex_state = 5},
  [153] = {.lex_state = 0},
  [154] = {.lex_state = 0},
  [155] = {.lex_state = 0},
  [156] = {.lex_state = 5},
  [157] = {.lex_state = 5},
  [158] = {.lex_state = 5},
  [159] = {.lex_state = 5},
  [160] = {.lex_state = 5},
  [161] = {.lex_state = 5},
  [162] = {.lex_state = 5},
  [163] = {.lex_state = 0},
  [164] = {.lex_state = 0},
  [165] = {.lex_state = 0},
  [166] = {.lex_state = 0},
  [167] = {.lex_state = 0},
  [168] = {.lex_state = 0},
  [169] = {.lex_state = 6},
  [170] = {.lex_state = 6},
  [171] = {.lex_state = 6},
  [172] = {.lex_state = 6},
  [173] = {.lex_state = 5},
  [174] = {.lex_state = 6},
  [175] = {.lex_state = 5},
  [176] = {.lex_state = 14},
  [177] = {.lex_state = 6},
  [178] = {.lex_state = 6},
  [179] = {.lex_state = 0},
  [180] = {.lex_state = 6},
  [181] = {.lex_state = 6},
  [182] = {.lex_state = 6},
  [183] = {.lex_state = 6},
  [184] = {.lex_state = 6},
  [185] = {.lex_state = 6},
  [186] = {.lex_state = 6},
  [187] = {.lex_state = 6},
  [188] = {.lex_state = 6},
  [189] = {.lex_state = 6},
  [190] = {.lex_state = 6},
  [191] = {.lex_state = 6},
  [192] = {.lex_state = 6},
  [193] = {.lex_state = 6},
  [194] = {.lex_state = 6},
  [195] = {.lex_state = 6},
  [196] = {.lex_state = 6},
  [197] = {.lex_state = 6},
  [198] = {.lex_state = 6},
  [199] = {.lex_state = 6},
  [200] = {.lex_state = 6},
  [201] = {.lex_state = 6},
  [202] = {.lex_state = 6},
  [203] = {.lex_state = 6},
  [204] = {.lex_state = 6},
  [205] = {.lex_state = 6},
  [206] = {.lex_state = 6},
  [207] = {.lex_state = 6},
  [208] = {.lex_state = 6},
  [209] = {.lex_state = 6},
  [210] = {.lex_state = 0},
  [211] = {.lex_state = 6},
  [212] = {.lex_state = 6},
  [213] = {.lex_state = 6},
  [214] = {.lex_state = 6},
  [215] = {.lex_state = 0},
  [216] = {.lex_state = 6},
  [217] = {.lex_state = 6},
  [218] = {.lex_state = 6},
  [219] = {.lex_state = 6},
  [220] = {.lex_state = 6},
  [221] = {.lex_state = 0},
  [222] = {.lex_state = 6},
  [223] = {.lex_state = 5},
  [224] = {.lex_state = 6},
  [225] = {.lex_state = 6},
  [226] = {.lex_state = 0},
  [227] = {.lex_state = 6},
  [228] = {.lex_state = 6},
  [229] = {.lex_state = 0},
  [230] = {.lex_state = 6},
  [231] = {.lex_state = 6},
  [232] = {.lex_state = 0},
  [233] = {.lex_state = 0},
  [234] = {.lex_state = 0},
  [235] = {.lex_state = 0},
  [236] = {.lex_state = 0},
  [237] = {.lex_state = 0},
  [238] = {.lex_state = 0},
  [239] = {.lex_state = 0},
  [240] = {.lex_state = 0},
  [241] = {.lex_state = 0},
  [242] = {.lex_state = 0},
  [243] = {.lex_state = 0},
  [244] = {.lex_state = 0},
  [245] = {.lex_state = 0},
  [246] = {.lex_state = 0},
  [247] = {.lex_state = 5},
  [248] = {.lex_state = 0},
  [249] = {.lex_state = 5},
  [250] = {.lex_state = 0},
  [251] = {.lex_state = 0},
  [252] = {.lex_state = 0},
  [253] = {.lex_state = 0},
  [254] = {.lex_state = 0},
  [255] = {.lex_state = 5},
  [256] = {.lex_state = 0},
  [257] = {.lex_state = 0},
  [258] = {.lex_state = 5},
  [259] = {.lex_state = 5},
  [260] = {.lex_state = 5},
  [261] = {.lex_state = 0},
  [262] = {.lex_state = 0},
  [263] = {.lex_state = 5},
  [264] = {.lex_state = 0},
  [265] = {.lex_state = 0},
  [266] = {.lex_state = 0},
  [267] = {.lex_state = 0},
  [268] = {.lex_state = 0},
  [269] = {.lex_state = 0},
  [270] = {.lex_state = 0},
  [271] = {.lex_state = 0},
  [272] = {.lex_state = 254},
  [273] = {.lex_state = 254},
  [274] = {.lex_state = 0},
  [275] = {.lex_state = 0},
  [276] = {.lex_state = 0},
  [277] = {.lex_state = 5},
  [278] = {.lex_state = 5},
  [279] = {.lex_state = 0},
  [280] = {.lex_state = 0},
  [281] = {.lex_state = 5},
  [282] = {.lex_state = 5},
  [283] = {.lex_state = 0},
  [284] = {.lex_state = 0},
  [285] = {.lex_state = 0},
  [286] = {.lex_state = 0},
  [287] = {.lex_state = 0},
  [288] = {.lex_state = 0},
  [289] = {.lex_state = 5},
  [290] = {.lex_state = 0},
  [291] = {.lex_state = 0},
  [292] = {.lex_state = 5},
  [293] = {.lex_state = 5},
  [294] = {.lex_state = 5},
  [295] = {.lex_state = 254},
  [296] = {.lex_state = 0},
  [297] = {.lex_state = 0},
  [298] = {.lex_state = 0},
  [299] = {.lex_state = 5},
  [300] = {.lex_state = 5},
  [301] = {.lex_state = 5},
  [302] = {.lex_state = 0},
  [303] = {.lex_state = 0},
  [304] = {.lex_state = 0},
  [305] = {.lex_state = 0},
  [306] = {.lex_state = 0},
  [307] = {.lex_state = 0},
  [308] = {.lex_state = 0},
  [309] = {.lex_state = 0},
  [310] = {.lex_state = 0},
  [311] = {.lex_state = 0},
  [312] = {.lex_state = 0},
  [313] = {.lex_state = 5},
  [314] = {.lex_state = 0},
  [315] = {.lex_state = 0},
  [316] = {.lex_state = 0},
  [317] = {.lex_state = 0},
  [318] = {.lex_state = 0},
  [319] = {.lex_state = 5},
  [320] = {.lex_state = 0},
  [321] = {.lex_state = 0},
  [322] = {.lex_state = 0},
  [323] = {.lex_state = 0},
  [324] = {.lex_state = 0},
  [325] = {.lex_state = 0},
  [326] = {.lex_state = 5},
  [327] = {.lex_state = 5},
  [328] = {.lex_state = 254},
  [329] = {.lex_state = 0},
  [330] = {.lex_state = 0},
  [331] = {.lex_state = 0},
  [332] = {.lex_state = 5},
  [333] = {.lex_state = 0},
  [334] = {.lex_state = 0},
  [335] = {.lex_state = 0},
  [336] = {.lex_state = 0},
  [337] = {.lex_state = 0},
};

static const uint16_t ts_parse_table[LARGE_STATE_COUNT][SYMBOL_COUNT] = {
  [0] = {
    [ts_builtin_sym_end] = ACTIONS(1),
    [anon_sym_POUNDinclude] = ACTIONS(1),
    [anon_sym_POUNDdefine] = ACTIONS(1),
    [anon_sym_POUNDifdef] = ACTIONS(1),
    [anon_sym_POUNDifndef] = ACTIONS(1),
    [anon_sym_POUNDendif] = ACTIONS(1),
    [anon_sym_EQ] = ACTIONS(1),
    [anon_sym_new] = ACTIONS(1),
    [anon_sym_for] = ACTIONS(1),
    [anon_sym_generate] = ACTIONS(1),
    [anon_sym_reaction] = ACTIONS(1),
    [anon_sym_collision] = ACTIONS(1),
    [anon_sym_LT_DASH_GT] = ACTIONS(1),
    [anon_sym_coulomb] = ACTIONS(1),
    [anon_sym_cross] = ACTIONS(1),
    [anon_sym_section] = ACTIONS(1),
    [anon_sym_metallic] = ACTIONS(1),
    [anon_sym_ks] = ACTIONS(1),
    [anon_sym_fermi_energy_ev] = ACTIONS(1),
    [anon_sym_ref_density] = ACTIONS(1),
    [anon_sym_excitation] = ACTIONS(1),
    [anon_sym_DASH_GT] = ACTIONS(1),
    [anon_sym_level] = ACTIONS(1),
    [anon_sym_get] = ACTIONS(1),
    [anon_sym_POUNDelse] = ACTIONS(1),
    [sym_define_key] = ACTIONS(1),
    [sym_define_ref] = ACTIONS(1),
    [sym_decimal] = ACTIONS(1),
    [anon_sym_true] = ACTIONS(1),
    [anon_sym_false] = ACTIONS(1),
    [anon_sym_yes] = ACTIONS(1),
    [anon_sym_no] = ACTIONS(1),
    [anon_sym_on] = ACTIONS(1),
    [anon_sym_off] = ACTIONS(1),
    [anon_sym_LBRACKdeg_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKrad_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKmrad_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKurad_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKum_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKmm_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKcm_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKm_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKfs_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKps_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKns_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKus_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKs_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACK_SLASHm3_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACK_SLASHcm3_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKkg_SLASHm3_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKg_SLASHcm3_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKJ_SLASHm3_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKJ_SLASHcm3_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKeV_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKK_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKPa_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKdynes_SLASHcm2_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKbar_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKergs_SLASHg_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKJ_SLASHkg_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKcm2_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKm2_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKcm2_SLASHs_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKm2_SLASHs_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKV_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKwebers_SLASHm_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKG_STARcm_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKV_SLASHm_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKV_SLASHcm_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKT_RBRACK] = ACTIONS(1),
    [anon_sym_LBRACKG_RBRACK] = ACTIONS(1),
    [anon_sym_1d] = ACTIONS(1),
    [anon_sym_2d] = ACTIONS(1),
    [anon_sym_3d] = ACTIONS(1),
    [anon_sym_SQUOTE] = ACTIONS(1),
    [anon_sym_DQUOTE] = ACTIONS(1),
    [anon_sym_LBRACE] = ACTIONS(1),
    [anon_sym_RBRACE] = ACTIONS(1),
    [anon_sym_LPAREN] = ACTIONS(1),
    [anon_sym_RPAREN] = ACTIONS(1),
    [sym_comment] = ACTIONS(3),
    [aux_sym__pterm_token1] = ACTIONS(1),
    [aux_sym__nterm_token1] = ACTIONS(1),
    [anon_sym_COLON] = ACTIONS(1),
    [anon_sym_rate] = ACTIONS(1),
    [anon_sym_janev_rate] = ACTIONS(1),
  },
  [1] = {
    [sym_input_file] = STATE(252),
    [sym__top] = STATE(18),
    [sym_include] = STATE(18),
    [sym_define] = STATE(18),
    [sym_ifxdef] = STATE(18),
    [sym__directive] = STATE(18),
    [sym_assignment] = STATE(18),
    [sym__statement] = STATE(18),
    [sym_new] = STATE(18),
    [sym_associative_new] = STATE(18),
    [sym_generate] = STATE(18),
    [sym_reaction] = STATE(18),
    [sym_collision] = STATE(18),
    [sym_excitation] = STATE(18),
    [sym_special_keys] = STATE(141),
    [sym_obj_key] = STATE(283),
    [aux_sym_input_file_repeat1] = STATE(18),
    [aux_sym_obj_key_repeat1] = STATE(141),
    [ts_builtin_sym_end] = ACTIONS(5),
    [anon_sym_POUNDinclude] = ACTIONS(7),
    [anon_sym_POUNDdefine] = ACTIONS(9),
    [anon_sym_POUNDifdef] = ACTIONS(11),
    [anon_sym_POUNDifndef] = ACTIONS(11),
    [anon_sym_new] = ACTIONS(13),
    [anon_sym_generate] = ACTIONS(15),
    [sym_identifier] = ACTIONS(17),
    [anon_sym_1d] = ACTIONS(19),
    [anon_sym_2d] = ACTIONS(19),
    [anon_sym_3d] = ACTIONS(19),
    [sym_comment] = ACTIONS(3),
  },
};

static const uint16_t ts_small_parse_table[] = {
  [0] = 5,
    ACTIONS(25), 1,
      aux_sym__define_value_token1,
    ACTIONS(27), 1,
      sym_comment,
    STATE(54), 1,
      sym_unit,
    ACTIONS(21), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
    ACTIONS(23), 37,
      anon_sym_LBRACKdeg_RBRACK,
      anon_sym_LBRACKrad_RBRACK,
      anon_sym_LBRACKmrad_RBRACK,
      anon_sym_LBRACKurad_RBRACK,
      anon_sym_LBRACKum_RBRACK,
      anon_sym_LBRACKmm_RBRACK,
      anon_sym_LBRACKcm_RBRACK,
      anon_sym_LBRACKm_RBRACK,
      anon_sym_LBRACKfs_RBRACK,
      anon_sym_LBRACKps_RBRACK,
      anon_sym_LBRACKns_RBRACK,
      anon_sym_LBRACKus_RBRACK,
      anon_sym_LBRACKs_RBRACK,
      anon_sym_LBRACK_SLASHm3_RBRACK,
      anon_sym_LBRACK_SLASHcm3_RBRACK,
      anon_sym_LBRACKkg_SLASHm3_RBRACK,
      anon_sym_LBRACKg_SLASHcm3_RBRACK,
      anon_sym_LBRACKJ_SLASHm3_RBRACK,
      anon_sym_LBRACKJ_SLASHcm3_RBRACK,
      anon_sym_LBRACKeV_RBRACK,
      anon_sym_LBRACKK_RBRACK,
      anon_sym_LBRACKPa_RBRACK,
      anon_sym_LBRACKdynes_SLASHcm2_RBRACK,
      anon_sym_LBRACKbar_RBRACK,
      anon_sym_LBRACKergs_SLASHg_RBRACK,
      anon_sym_LBRACKJ_SLASHkg_RBRACK,
      anon_sym_LBRACKcm2_RBRACK,
      anon_sym_LBRACKm2_RBRACK,
      anon_sym_LBRACKcm2_SLASHs_RBRACK,
      anon_sym_LBRACKm2_SLASHs_RBRACK,
      anon_sym_LBRACKV_RBRACK,
      anon_sym_LBRACKwebers_SLASHm_RBRACK,
      anon_sym_LBRACKG_STARcm_RBRACK,
      anon_sym_LBRACKV_SLASHm_RBRACK,
      anon_sym_LBRACKV_SLASHcm_RBRACK,
      anon_sym_LBRACKT_RBRACK,
      anon_sym_LBRACKG_RBRACK,
  [69] = 5,
    ACTIONS(3), 1,
      sym_comment,
    STATE(109), 1,
      sym_unit,
    ACTIONS(25), 6,
      sym_define_ref,
      sym_decimal,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_RBRACE,
      anon_sym_RPAREN,
    ACTIONS(21), 7,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
    ACTIONS(29), 37,
      anon_sym_LBRACKdeg_RBRACK,
      anon_sym_LBRACKrad_RBRACK,
      anon_sym_LBRACKmrad_RBRACK,
      anon_sym_LBRACKurad_RBRACK,
      anon_sym_LBRACKum_RBRACK,
      anon_sym_LBRACKmm_RBRACK,
      anon_sym_LBRACKcm_RBRACK,
      anon_sym_LBRACKm_RBRACK,
      anon_sym_LBRACKfs_RBRACK,
      anon_sym_LBRACKps_RBRACK,
      anon_sym_LBRACKns_RBRACK,
      anon_sym_LBRACKus_RBRACK,
      anon_sym_LBRACKs_RBRACK,
      anon_sym_LBRACK_SLASHm3_RBRACK,
      anon_sym_LBRACK_SLASHcm3_RBRACK,
      anon_sym_LBRACKkg_SLASHm3_RBRACK,
      anon_sym_LBRACKg_SLASHcm3_RBRACK,
      anon_sym_LBRACKJ_SLASHm3_RBRACK,
      anon_sym_LBRACKJ_SLASHcm3_RBRACK,
      anon_sym_LBRACKeV_RBRACK,
      anon_sym_LBRACKK_RBRACK,
      anon_sym_LBRACKPa_RBRACK,
      anon_sym_LBRACKdynes_SLASHcm2_RBRACK,
      anon_sym_LBRACKbar_RBRACK,
      anon_sym_LBRACKergs_SLASHg_RBRACK,
      anon_sym_LBRACKJ_SLASHkg_RBRACK,
      anon_sym_LBRACKcm2_RBRACK,
      anon_sym_LBRACKm2_RBRACK,
      anon_sym_LBRACKcm2_SLASHs_RBRACK,
      anon_sym_LBRACKm2_SLASHs_RBRACK,
      anon_sym_LBRACKV_RBRACK,
      anon_sym_LBRACKwebers_SLASHm_RBRACK,
      anon_sym_LBRACKG_STARcm_RBRACK,
      anon_sym_LBRACKV_SLASHm_RBRACK,
      anon_sym_LBRACKV_SLASHcm_RBRACK,
      anon_sym_LBRACKT_RBRACK,
      anon_sym_LBRACKG_RBRACK,
  [132] = 5,
    ACTIONS(3), 1,
      sym_comment,
    STATE(85), 1,
      sym_unit,
    ACTIONS(21), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(25), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(31), 37,
      anon_sym_LBRACKdeg_RBRACK,
      anon_sym_LBRACKrad_RBRACK,
      anon_sym_LBRACKmrad_RBRACK,
      anon_sym_LBRACKurad_RBRACK,
      anon_sym_LBRACKum_RBRACK,
      anon_sym_LBRACKmm_RBRACK,
      anon_sym_LBRACKcm_RBRACK,
      anon_sym_LBRACKm_RBRACK,
      anon_sym_LBRACKfs_RBRACK,
      anon_sym_LBRACKps_RBRACK,
      anon_sym_LBRACKns_RBRACK,
      anon_sym_LBRACKus_RBRACK,
      anon_sym_LBRACKs_RBRACK,
      anon_sym_LBRACK_SLASHm3_RBRACK,
      anon_sym_LBRACK_SLASHcm3_RBRACK,
      anon_sym_LBRACKkg_SLASHm3_RBRACK,
      anon_sym_LBRACKg_SLASHcm3_RBRACK,
      anon_sym_LBRACKJ_SLASHm3_RBRACK,
      anon_sym_LBRACKJ_SLASHcm3_RBRACK,
      anon_sym_LBRACKeV_RBRACK,
      anon_sym_LBRACKK_RBRACK,
      anon_sym_LBRACKPa_RBRACK,
      anon_sym_LBRACKdynes_SLASHcm2_RBRACK,
      anon_sym_LBRACKbar_RBRACK,
      anon_sym_LBRACKergs_SLASHg_RBRACK,
      anon_sym_LBRACKJ_SLASHkg_RBRACK,
      anon_sym_LBRACKcm2_RBRACK,
      anon_sym_LBRACKm2_RBRACK,
      anon_sym_LBRACKcm2_SLASHs_RBRACK,
      anon_sym_LBRACKm2_SLASHs_RBRACK,
      anon_sym_LBRACKV_RBRACK,
      anon_sym_LBRACKwebers_SLASHm_RBRACK,
      anon_sym_LBRACKG_STARcm_RBRACK,
      anon_sym_LBRACKV_SLASHm_RBRACK,
      anon_sym_LBRACKV_SLASHcm_RBRACK,
      anon_sym_LBRACKT_RBRACK,
      anon_sym_LBRACKG_RBRACK,
  [195] = 5,
    ACTIONS(3), 1,
      sym_comment,
    STATE(118), 1,
      sym_unit,
    ACTIONS(21), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(25), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
    ACTIONS(33), 37,
      anon_sym_LBRACKdeg_RBRACK,
      anon_sym_LBRACKrad_RBRACK,
      anon_sym_LBRACKmrad_RBRACK,
      anon_sym_LBRACKurad_RBRACK,
      anon_sym_LBRACKum_RBRACK,
      anon_sym_LBRACKmm_RBRACK,
      anon_sym_LBRACKcm_RBRACK,
      anon_sym_LBRACKm_RBRACK,
      anon_sym_LBRACKfs_RBRACK,
      anon_sym_LBRACKps_RBRACK,
      anon_sym_LBRACKns_RBRACK,
      anon_sym_LBRACKus_RBRACK,
      anon_sym_LBRACKs_RBRACK,
      anon_sym_LBRACK_SLASHm3_RBRACK,
      anon_sym_LBRACK_SLASHcm3_RBRACK,
      anon_sym_LBRACKkg_SLASHm3_RBRACK,
      anon_sym_LBRACKg_SLASHcm3_RBRACK,
      anon_sym_LBRACKJ_SLASHm3_RBRACK,
      anon_sym_LBRACKJ_SLASHcm3_RBRACK,
      anon_sym_LBRACKeV_RBRACK,
      anon_sym_LBRACKK_RBRACK,
      anon_sym_LBRACKPa_RBRACK,
      anon_sym_LBRACKdynes_SLASHcm2_RBRACK,
      anon_sym_LBRACKbar_RBRACK,
      anon_sym_LBRACKergs_SLASHg_RBRACK,
      anon_sym_LBRACKJ_SLASHkg_RBRACK,
      anon_sym_LBRACKcm2_RBRACK,
      anon_sym_LBRACKm2_RBRACK,
      anon_sym_LBRACKcm2_SLASHs_RBRACK,
      anon_sym_LBRACKm2_SLASHs_RBRACK,
      anon_sym_LBRACKV_RBRACK,
      anon_sym_LBRACKwebers_SLASHm_RBRACK,
      anon_sym_LBRACKG_STARcm_RBRACK,
      anon_sym_LBRACKV_SLASHm_RBRACK,
      anon_sym_LBRACKV_SLASHcm_RBRACK,
      anon_sym_LBRACKT_RBRACK,
      anon_sym_LBRACKG_RBRACK,
  [256] = 18,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(35), 1,
      anon_sym_new,
    ACTIONS(37), 1,
      anon_sym_generate,
    ACTIONS(39), 1,
      sym_define_ref,
    ACTIONS(41), 1,
      sym_decimal,
    ACTIONS(45), 1,
      sym_identifier,
    ACTIONS(49), 1,
      anon_sym_SQUOTE,
    ACTIONS(51), 1,
      anon_sym_DQUOTE,
    ACTIONS(53), 1,
      anon_sym_LBRACE,
    ACTIONS(55), 1,
      anon_sym_LPAREN,
    ACTIONS(57), 1,
      aux_sym__define_value_token1,
    STATE(34), 1,
      sym_obj_key,
    STATE(133), 1,
      sym__define_value,
    STATE(28), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    STATE(43), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    ACTIONS(47), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(43), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
    STATE(10), 17,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      sym_block,
      sym_tuple,
      sym_list,
      aux_sym__define_value_repeat1,
  [336] = 18,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(35), 1,
      anon_sym_new,
    ACTIONS(37), 1,
      anon_sym_generate,
    ACTIONS(41), 1,
      sym_decimal,
    ACTIONS(45), 1,
      sym_identifier,
    ACTIONS(49), 1,
      anon_sym_SQUOTE,
    ACTIONS(51), 1,
      anon_sym_DQUOTE,
    ACTIONS(53), 1,
      anon_sym_LBRACE,
    ACTIONS(55), 1,
      anon_sym_LPAREN,
    ACTIONS(59), 1,
      sym_define_ref,
    ACTIONS(61), 1,
      aux_sym__define_value_token1,
    STATE(34), 1,
      sym_obj_key,
    STATE(80), 1,
      sym__define_value,
    STATE(28), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    STATE(43), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    ACTIONS(47), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(43), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
    STATE(8), 17,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      sym_block,
      sym_tuple,
      sym_list,
      aux_sym__define_value_repeat1,
  [416] = 17,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(35), 1,
      anon_sym_new,
    ACTIONS(37), 1,
      anon_sym_generate,
    ACTIONS(41), 1,
      sym_decimal,
    ACTIONS(45), 1,
      sym_identifier,
    ACTIONS(49), 1,
      anon_sym_SQUOTE,
    ACTIONS(51), 1,
      anon_sym_DQUOTE,
    ACTIONS(53), 1,
      anon_sym_LBRACE,
    ACTIONS(55), 1,
      anon_sym_LPAREN,
    ACTIONS(63), 1,
      sym_define_ref,
    ACTIONS(65), 1,
      aux_sym__define_value_token1,
    STATE(34), 1,
      sym_obj_key,
    STATE(28), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    STATE(43), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    ACTIONS(47), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(43), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
    STATE(11), 17,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      sym_block,
      sym_tuple,
      sym_list,
      aux_sym__define_value_repeat1,
  [493] = 19,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(67), 1,
      anon_sym_POUNDinclude,
    ACTIONS(69), 1,
      anon_sym_POUNDdefine,
    ACTIONS(73), 1,
      anon_sym_new,
    ACTIONS(75), 1,
      anon_sym_get,
    ACTIONS(77), 1,
      sym_define_ref,
    ACTIONS(79), 1,
      sym_decimal,
    ACTIONS(83), 1,
      sym_identifier,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(89), 1,
      anon_sym_RBRACE,
    STATE(276), 1,
      sym_obj_key,
    ACTIONS(71), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(59), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
    STATE(20), 10,
      sym__nested,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__nested_directive,
      sym_assignment,
      sym__nested_statement,
      sym_new,
      sym_get,
      aux_sym_block_repeat1,
  [574] = 17,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(35), 1,
      anon_sym_new,
    ACTIONS(37), 1,
      anon_sym_generate,
    ACTIONS(41), 1,
      sym_decimal,
    ACTIONS(45), 1,
      sym_identifier,
    ACTIONS(49), 1,
      anon_sym_SQUOTE,
    ACTIONS(51), 1,
      anon_sym_DQUOTE,
    ACTIONS(53), 1,
      anon_sym_LBRACE,
    ACTIONS(55), 1,
      anon_sym_LPAREN,
    ACTIONS(63), 1,
      sym_define_ref,
    ACTIONS(91), 1,
      aux_sym__define_value_token1,
    STATE(34), 1,
      sym_obj_key,
    STATE(28), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    STATE(43), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    ACTIONS(47), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(43), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
    STATE(11), 17,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      sym_block,
      sym_tuple,
      sym_list,
      aux_sym__define_value_repeat1,
  [651] = 17,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(93), 1,
      anon_sym_new,
    ACTIONS(96), 1,
      anon_sym_generate,
    ACTIONS(99), 1,
      sym_define_ref,
    ACTIONS(102), 1,
      sym_decimal,
    ACTIONS(108), 1,
      sym_identifier,
    ACTIONS(114), 1,
      anon_sym_SQUOTE,
    ACTIONS(117), 1,
      anon_sym_DQUOTE,
    ACTIONS(120), 1,
      anon_sym_LBRACE,
    ACTIONS(123), 1,
      anon_sym_LPAREN,
    ACTIONS(126), 1,
      aux_sym__define_value_token1,
    STATE(34), 1,
      sym_obj_key,
    STATE(28), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    STATE(43), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    ACTIONS(111), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(105), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
    STATE(11), 17,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      sym_block,
      sym_tuple,
      sym_list,
      aux_sym__define_value_repeat1,
  [728] = 14,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(7), 1,
      anon_sym_POUNDinclude,
    ACTIONS(9), 1,
      anon_sym_POUNDdefine,
    ACTIONS(13), 1,
      anon_sym_new,
    ACTIONS(15), 1,
      anon_sym_generate,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(128), 1,
      anon_sym_POUNDendif,
    ACTIONS(130), 1,
      anon_sym_POUNDelse,
    STATE(254), 1,
      sym_else_block,
    STATE(283), 1,
      sym_obj_key,
    ACTIONS(11), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(15), 14,
      sym__top,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      aux_sym_input_file_repeat1,
  [788] = 14,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(7), 1,
      anon_sym_POUNDinclude,
    ACTIONS(9), 1,
      anon_sym_POUNDdefine,
    ACTIONS(13), 1,
      anon_sym_new,
    ACTIONS(15), 1,
      anon_sym_generate,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(130), 1,
      anon_sym_POUNDelse,
    ACTIONS(132), 1,
      anon_sym_POUNDendif,
    STATE(283), 1,
      sym_obj_key,
    STATE(302), 1,
      sym_else_block,
    ACTIONS(11), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(12), 14,
      sym__top,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      aux_sym_input_file_repeat1,
  [848] = 14,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(7), 1,
      anon_sym_POUNDinclude,
    ACTIONS(9), 1,
      anon_sym_POUNDdefine,
    ACTIONS(13), 1,
      anon_sym_new,
    ACTIONS(15), 1,
      anon_sym_generate,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(130), 1,
      anon_sym_POUNDelse,
    ACTIONS(134), 1,
      anon_sym_POUNDendif,
    STATE(253), 1,
      sym_else_block,
    STATE(283), 1,
      sym_obj_key,
    ACTIONS(11), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(15), 14,
      sym__top,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      aux_sym_input_file_repeat1,
  [908] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(138), 1,
      anon_sym_POUNDinclude,
    ACTIONS(141), 1,
      anon_sym_POUNDdefine,
    ACTIONS(147), 1,
      anon_sym_new,
    ACTIONS(150), 1,
      anon_sym_generate,
    ACTIONS(153), 1,
      sym_identifier,
    STATE(283), 1,
      sym_obj_key,
    ACTIONS(144), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(136), 3,
      ts_builtin_sym_end,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
    ACTIONS(156), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(15), 14,
      sym__top,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      aux_sym_input_file_repeat1,
  [964] = 14,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(7), 1,
      anon_sym_POUNDinclude,
    ACTIONS(9), 1,
      anon_sym_POUNDdefine,
    ACTIONS(13), 1,
      anon_sym_new,
    ACTIONS(15), 1,
      anon_sym_generate,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(130), 1,
      anon_sym_POUNDelse,
    ACTIONS(159), 1,
      anon_sym_POUNDendif,
    STATE(264), 1,
      sym_else_block,
    STATE(283), 1,
      sym_obj_key,
    ACTIONS(11), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(14), 14,
      sym__top,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      aux_sym_input_file_repeat1,
  [1024] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(7), 1,
      anon_sym_POUNDinclude,
    ACTIONS(9), 1,
      anon_sym_POUNDdefine,
    ACTIONS(13), 1,
      anon_sym_new,
    ACTIONS(15), 1,
      anon_sym_generate,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(161), 1,
      anon_sym_POUNDendif,
    STATE(283), 1,
      sym_obj_key,
    ACTIONS(11), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(15), 14,
      sym__top,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      aux_sym_input_file_repeat1,
  [1078] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(7), 1,
      anon_sym_POUNDinclude,
    ACTIONS(9), 1,
      anon_sym_POUNDdefine,
    ACTIONS(13), 1,
      anon_sym_new,
    ACTIONS(15), 1,
      anon_sym_generate,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(163), 1,
      ts_builtin_sym_end,
    STATE(283), 1,
      sym_obj_key,
    ACTIONS(11), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(15), 14,
      sym__top,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      aux_sym_input_file_repeat1,
  [1132] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(7), 1,
      anon_sym_POUNDinclude,
    ACTIONS(9), 1,
      anon_sym_POUNDdefine,
    ACTIONS(13), 1,
      anon_sym_new,
    ACTIONS(15), 1,
      anon_sym_generate,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(165), 1,
      anon_sym_POUNDendif,
    STATE(283), 1,
      sym_obj_key,
    ACTIONS(11), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(17), 14,
      sym__top,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__directive,
      sym_assignment,
      sym__statement,
      sym_new,
      sym_associative_new,
      sym_generate,
      sym_reaction,
      sym_collision,
      sym_excitation,
      aux_sym_input_file_repeat1,
  [1186] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(67), 1,
      anon_sym_POUNDinclude,
    ACTIONS(69), 1,
      anon_sym_POUNDdefine,
    ACTIONS(73), 1,
      anon_sym_new,
    ACTIONS(75), 1,
      anon_sym_get,
    ACTIONS(167), 1,
      anon_sym_RBRACE,
    STATE(276), 1,
      sym_obj_key,
    ACTIONS(71), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(24), 10,
      sym__nested,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__nested_directive,
      sym_assignment,
      sym__nested_statement,
      sym_new,
      sym_get,
      aux_sym_block_repeat1,
  [1236] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(67), 1,
      anon_sym_POUNDinclude,
    ACTIONS(69), 1,
      anon_sym_POUNDdefine,
    ACTIONS(73), 1,
      anon_sym_new,
    ACTIONS(75), 1,
      anon_sym_get,
    ACTIONS(169), 1,
      anon_sym_RBRACE,
    STATE(276), 1,
      sym_obj_key,
    ACTIONS(71), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(24), 10,
      sym__nested,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__nested_directive,
      sym_assignment,
      sym__nested_statement,
      sym_new,
      sym_get,
      aux_sym_block_repeat1,
  [1286] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(67), 1,
      anon_sym_POUNDinclude,
    ACTIONS(69), 1,
      anon_sym_POUNDdefine,
    ACTIONS(73), 1,
      anon_sym_new,
    ACTIONS(75), 1,
      anon_sym_get,
    ACTIONS(171), 1,
      anon_sym_RBRACE,
    STATE(276), 1,
      sym_obj_key,
    ACTIONS(71), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(25), 10,
      sym__nested,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__nested_directive,
      sym_assignment,
      sym__nested_statement,
      sym_new,
      sym_get,
      aux_sym_block_repeat1,
  [1336] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(67), 1,
      anon_sym_POUNDinclude,
    ACTIONS(69), 1,
      anon_sym_POUNDdefine,
    ACTIONS(73), 1,
      anon_sym_new,
    ACTIONS(75), 1,
      anon_sym_get,
    ACTIONS(173), 1,
      anon_sym_RBRACE,
    STATE(276), 1,
      sym_obj_key,
    ACTIONS(71), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(21), 10,
      sym__nested,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__nested_directive,
      sym_assignment,
      sym__nested_statement,
      sym_new,
      sym_get,
      aux_sym_block_repeat1,
  [1386] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(175), 1,
      anon_sym_POUNDinclude,
    ACTIONS(178), 1,
      anon_sym_POUNDdefine,
    ACTIONS(184), 1,
      anon_sym_new,
    ACTIONS(187), 1,
      anon_sym_get,
    ACTIONS(190), 1,
      sym_identifier,
    ACTIONS(196), 1,
      anon_sym_RBRACE,
    STATE(276), 1,
      sym_obj_key,
    ACTIONS(181), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(193), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(24), 10,
      sym__nested,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__nested_directive,
      sym_assignment,
      sym__nested_statement,
      sym_new,
      sym_get,
      aux_sym_block_repeat1,
  [1436] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(67), 1,
      anon_sym_POUNDinclude,
    ACTIONS(69), 1,
      anon_sym_POUNDdefine,
    ACTIONS(73), 1,
      anon_sym_new,
    ACTIONS(75), 1,
      anon_sym_get,
    ACTIONS(198), 1,
      anon_sym_RBRACE,
    STATE(276), 1,
      sym_obj_key,
    ACTIONS(71), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(24), 10,
      sym__nested,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__nested_directive,
      sym_assignment,
      sym__nested_statement,
      sym_new,
      sym_get,
      aux_sym_block_repeat1,
  [1486] = 12,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(17), 1,
      sym_identifier,
    ACTIONS(67), 1,
      anon_sym_POUNDinclude,
    ACTIONS(69), 1,
      anon_sym_POUNDdefine,
    ACTIONS(73), 1,
      anon_sym_new,
    ACTIONS(75), 1,
      anon_sym_get,
    ACTIONS(89), 1,
      anon_sym_RBRACE,
    STATE(276), 1,
      sym_obj_key,
    ACTIONS(71), 2,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    STATE(20), 10,
      sym__nested,
      sym_include,
      sym_define,
      sym_ifxdef,
      sym__nested_directive,
      sym_assignment,
      sym__nested_statement,
      sym_new,
      sym_get,
      aux_sym_block_repeat1,
  [1536] = 6,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(202), 1,
      sym_identifier,
    ACTIONS(208), 1,
      aux_sym__define_value_token1,
    STATE(27), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(205), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(200), 15,
      anon_sym_EQ,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [1572] = 4,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(212), 1,
      aux_sym__define_value_token1,
    STATE(27), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(210), 19,
      anon_sym_EQ,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [1604] = 11,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(214), 1,
      sym_define_ref,
    ACTIONS(216), 1,
      sym_decimal,
    ACTIONS(218), 1,
      sym_identifier,
    ACTIONS(220), 1,
      anon_sym_SQUOTE,
    ACTIONS(222), 1,
      anon_sym_DQUOTE,
    ACTIONS(224), 1,
      anon_sym_LBRACE,
    ACTIONS(226), 1,
      anon_sym_LPAREN,
    STATE(43), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    ACTIONS(43), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
    STATE(56), 6,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      sym_tuple,
      sym_list,
  [1649] = 11,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(228), 1,
      sym_define_ref,
    ACTIONS(230), 1,
      sym_decimal,
    ACTIONS(234), 1,
      sym_identifier,
    ACTIONS(236), 1,
      anon_sym_SQUOTE,
    ACTIONS(238), 1,
      anon_sym_DQUOTE,
    ACTIONS(240), 1,
      anon_sym_LBRACE,
    ACTIONS(242), 1,
      anon_sym_LPAREN,
    STATE(77), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    ACTIONS(232), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
    STATE(116), 6,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      sym_tuple,
      sym_list,
  [1694] = 11,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(244), 1,
      sym_define_ref,
    ACTIONS(246), 1,
      sym_decimal,
    ACTIONS(250), 1,
      sym_identifier,
    ACTIONS(252), 1,
      anon_sym_SQUOTE,
    ACTIONS(254), 1,
      anon_sym_DQUOTE,
    ACTIONS(256), 1,
      anon_sym_LBRACE,
    ACTIONS(258), 1,
      anon_sym_LPAREN,
    STATE(134), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    ACTIONS(248), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
    STATE(120), 6,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      sym_tuple,
      sym_list,
  [1739] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(262), 1,
      aux_sym__define_value_token1,
    ACTIONS(260), 19,
      anon_sym_EQ,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [1767] = 4,
    ACTIONS(25), 1,
      aux_sym__define_value_token1,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(264), 1,
      anon_sym_EQ,
    ACTIONS(21), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [1797] = 4,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(266), 1,
      anon_sym_EQ,
    ACTIONS(270), 1,
      aux_sym__define_value_token1,
    ACTIONS(268), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [1827] = 10,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(272), 1,
      sym_define_ref,
    ACTIONS(275), 1,
      sym_decimal,
    ACTIONS(281), 1,
      sym_identifier,
    ACTIONS(284), 1,
      anon_sym_SQUOTE,
    ACTIONS(287), 1,
      anon_sym_DQUOTE,
    ACTIONS(290), 2,
      anon_sym_RBRACE,
      anon_sym_RPAREN,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(35), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(278), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [1869] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(294), 1,
      aux_sym__define_value_token1,
    ACTIONS(292), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [1896] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(298), 1,
      aux_sym__define_value_token1,
    ACTIONS(296), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [1923] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(302), 1,
      aux_sym__define_value_token1,
    ACTIONS(300), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [1950] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(306), 1,
      aux_sym__define_value_token1,
    ACTIONS(304), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [1977] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(310), 1,
      aux_sym__define_value_token1,
    ACTIONS(308), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2004] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(314), 1,
      aux_sym__define_value_token1,
    ACTIONS(312), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2031] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(318), 1,
      aux_sym__define_value_token1,
    ACTIONS(316), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2058] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(322), 1,
      aux_sym__define_value_token1,
    ACTIONS(320), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2085] = 10,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(324), 1,
      sym_define_ref,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(328), 1,
      sym_identifier,
    ACTIONS(330), 1,
      anon_sym_RBRACE,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(35), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [2126] = 10,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(324), 1,
      sym_define_ref,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(328), 1,
      sym_identifier,
    ACTIONS(332), 1,
      anon_sym_RPAREN,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(35), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [2167] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(336), 1,
      aux_sym__define_value_token1,
    ACTIONS(334), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2194] = 10,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(324), 1,
      sym_define_ref,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(328), 1,
      sym_identifier,
    ACTIONS(338), 1,
      anon_sym_RPAREN,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(35), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [2235] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(342), 1,
      aux_sym__define_value_token1,
    ACTIONS(340), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2262] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(346), 1,
      aux_sym__define_value_token1,
    ACTIONS(344), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2289] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(350), 1,
      aux_sym__define_value_token1,
    ACTIONS(348), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2316] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(354), 1,
      aux_sym__define_value_token1,
    ACTIONS(352), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2343] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(358), 1,
      aux_sym__define_value_token1,
    ACTIONS(356), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2370] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(362), 1,
      aux_sym__define_value_token1,
    ACTIONS(360), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2397] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(366), 1,
      aux_sym__define_value_token1,
    ACTIONS(364), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2424] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(370), 1,
      aux_sym__define_value_token1,
    ACTIONS(368), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2451] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(374), 1,
      aux_sym__define_value_token1,
    ACTIONS(372), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2478] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(378), 1,
      aux_sym__define_value_token1,
    ACTIONS(376), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2505] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(382), 1,
      aux_sym__define_value_token1,
    ACTIONS(380), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2532] = 10,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(324), 1,
      sym_define_ref,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(328), 1,
      sym_identifier,
    ACTIONS(384), 1,
      anon_sym_RBRACE,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(35), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [2573] = 10,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(324), 1,
      sym_define_ref,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(328), 1,
      sym_identifier,
    ACTIONS(386), 1,
      anon_sym_RPAREN,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(35), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [2614] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(390), 1,
      aux_sym__define_value_token1,
    ACTIONS(388), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2641] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(394), 1,
      aux_sym__define_value_token1,
    ACTIONS(392), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2668] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(398), 1,
      aux_sym__define_value_token1,
    ACTIONS(396), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2695] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(402), 1,
      aux_sym__define_value_token1,
    ACTIONS(400), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2722] = 10,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(324), 1,
      sym_define_ref,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(328), 1,
      sym_identifier,
    ACTIONS(404), 1,
      anon_sym_RBRACE,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(35), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [2763] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(408), 1,
      aux_sym__define_value_token1,
    ACTIONS(406), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2790] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(412), 1,
      aux_sym__define_value_token1,
    ACTIONS(410), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2817] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(416), 1,
      aux_sym__define_value_token1,
    ACTIONS(414), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2844] = 3,
    ACTIONS(27), 1,
      sym_comment,
    ACTIONS(420), 1,
      aux_sym__define_value_token1,
    ACTIONS(418), 18,
      anon_sym_new,
      anon_sym_generate,
      sym_define_ref,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
      anon_sym_LPAREN,
  [2871] = 9,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(77), 1,
      sym_define_ref,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(422), 1,
      sym_identifier,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(59), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [2909] = 9,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(424), 1,
      sym_define_ref,
    ACTIONS(426), 1,
      sym_identifier,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(47), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [2947] = 9,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(428), 1,
      sym_define_ref,
    ACTIONS(430), 1,
      sym_identifier,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(65), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [2985] = 9,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(432), 1,
      sym_define_ref,
    ACTIONS(434), 1,
      sym_identifier,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(45), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [3023] = 9,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(436), 1,
      sym_define_ref,
    ACTIONS(438), 1,
      sym_identifier,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(60), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [3061] = 9,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(85), 1,
      anon_sym_SQUOTE,
    ACTIONS(87), 1,
      anon_sym_DQUOTE,
    ACTIONS(326), 1,
      sym_decimal,
    ACTIONS(440), 1,
      sym_define_ref,
    ACTIONS(442), 1,
      sym_identifier,
    STATE(98), 2,
      sym__string_literal_single,
      sym__string_literal_double,
    STATE(44), 5,
      sym_boolean,
      sym_dimension,
      sym_string_literal,
      sym__value,
      aux_sym_tuple_repeat1,
    ACTIONS(81), 6,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
  [3099] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(25), 4,
      sym_define_ref,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_RBRACE,
    ACTIONS(444), 4,
      anon_sym_EQ,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(21), 8,
      sym_decimal,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
  [3125] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(320), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(322), 11,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_LBRACE,
  [3147] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(348), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(350), 11,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_LBRACE,
  [3169] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(312), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(314), 11,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_LBRACE,
  [3191] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(448), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(446), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3212] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(400), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(402), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3233] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(296), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(298), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3254] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(356), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(358), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3275] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(344), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(346), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3296] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(364), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(366), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3317] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(316), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(318), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3338] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(406), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(408), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3359] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(300), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(302), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3380] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(308), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(310), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3401] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(418), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(420), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3422] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(304), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(306), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3443] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(388), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(390), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3464] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(340), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(342), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3485] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(350), 6,
      sym_define_ref,
      sym_decimal,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_RBRACE,
      anon_sym_RPAREN,
    ACTIONS(348), 7,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
  [3506] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(314), 6,
      sym_define_ref,
      sym_decimal,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_RBRACE,
      anon_sym_RPAREN,
    ACTIONS(312), 7,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
  [3527] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(452), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(450), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3548] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(352), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(354), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3569] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(322), 6,
      sym_define_ref,
      sym_decimal,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_RBRACE,
      anon_sym_RPAREN,
    ACTIONS(320), 7,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
  [3590] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(360), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(362), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3611] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(334), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(336), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3632] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(456), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(454), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3653] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(414), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(416), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3674] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(376), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(378), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3695] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(460), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(458), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3716] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(396), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(398), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3737] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(292), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(294), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3758] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(464), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(462), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3779] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(392), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(394), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3800] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(366), 6,
      sym_define_ref,
      sym_decimal,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_RBRACE,
      anon_sym_RPAREN,
    ACTIONS(364), 7,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
  [3821] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(362), 6,
      sym_define_ref,
      sym_decimal,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_RBRACE,
      anon_sym_RPAREN,
    ACTIONS(360), 7,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
  [3842] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(410), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(412), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3863] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(380), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(382), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3884] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(368), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(370), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3905] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(378), 6,
      sym_define_ref,
      sym_decimal,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_RBRACE,
      anon_sym_RPAREN,
    ACTIONS(376), 7,
      anon_sym_true,
      anon_sym_false,
      anon_sym_yes,
      anon_sym_no,
      anon_sym_on,
      anon_sym_off,
      sym_identifier,
  [3926] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(468), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(466), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3947] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(372), 3,
      anon_sym_new,
      anon_sym_generate,
      sym_identifier,
    ACTIONS(374), 10,
      ts_builtin_sym_end,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_POUNDendif,
      anon_sym_POUNDelse,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [3968] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(460), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(458), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [3987] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(364), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(366), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4006] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(376), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(378), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4025] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(372), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(374), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4044] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(296), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(298), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4063] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(472), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(470), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4082] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(352), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(354), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4101] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(456), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(454), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4120] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(468), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(466), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4139] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(356), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(358), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4158] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(344), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(346), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4177] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(300), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(302), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4196] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(392), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(394), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4215] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(348), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(350), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4234] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(464), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(462), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4253] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(360), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(362), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4272] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(448), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(446), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4291] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(320), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(322), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4310] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(452), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(450), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4329] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(312), 3,
      anon_sym_new,
      anon_sym_get,
      sym_identifier,
    ACTIONS(314), 8,
      anon_sym_POUNDinclude,
      anon_sym_POUNDdefine,
      anon_sym_POUNDifdef,
      anon_sym_POUNDifndef,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_RBRACE,
  [4348] = 6,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(210), 1,
      anon_sym_for,
    ACTIONS(474), 1,
      sym_identifier,
    STATE(140), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(212), 3,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
    ACTIONS(476), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [4372] = 8,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(478), 1,
      anon_sym_reaction,
    ACTIONS(480), 1,
      anon_sym_collision,
    ACTIONS(482), 1,
      anon_sym_excitation,
    ACTIONS(484), 1,
      sym_identifier,
    STATE(147), 1,
      sym_obj_key,
    STATE(137), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(476), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [4400] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(486), 1,
      sym_identifier,
    STATE(139), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(489), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(208), 4,
      anon_sym_EQ,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
  [4422] = 6,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(200), 1,
      anon_sym_for,
    ACTIONS(492), 1,
      sym_identifier,
    STATE(140), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(208), 3,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
    ACTIONS(495), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [4446] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(498), 1,
      sym_identifier,
    STATE(139), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
    ACTIONS(212), 4,
      anon_sym_EQ,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
  [4468] = 8,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(484), 1,
      sym_identifier,
    ACTIONS(500), 1,
      anon_sym_reaction,
    ACTIONS(502), 1,
      anon_sym_collision,
    ACTIONS(504), 1,
      anon_sym_excitation,
    STATE(143), 1,
      sym_obj_key,
    STATE(137), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(476), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [4496] = 8,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(506), 1,
      anon_sym_for,
    ACTIONS(508), 1,
      anon_sym_SQUOTE,
    ACTIONS(510), 1,
      anon_sym_DQUOTE,
    ACTIONS(512), 1,
      anon_sym_LBRACE,
    STATE(62), 1,
      sym_block,
    STATE(179), 1,
      sym_string_literal,
    STATE(234), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4522] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(262), 8,
      anon_sym_EQ,
      sym_identifier,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
  [4536] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(516), 1,
      aux_sym__pterm_token1,
    ACTIONS(519), 1,
      aux_sym__nterm_token1,
    ACTIONS(514), 3,
      anon_sym_DASH_GT,
      anon_sym_RBRACE,
      anon_sym_COLON,
    STATE(145), 3,
      sym__pterm,
      sym__nterm,
      aux_sym_chems_repeat1,
  [4556] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(524), 1,
      aux_sym__pterm_token1,
    ACTIONS(526), 1,
      aux_sym__nterm_token1,
    ACTIONS(522), 3,
      anon_sym_DASH_GT,
      anon_sym_RBRACE,
      anon_sym_COLON,
    STATE(145), 3,
      sym__pterm,
      sym__nterm,
      aux_sym_chems_repeat1,
  [4576] = 8,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(508), 1,
      anon_sym_SQUOTE,
    ACTIONS(510), 1,
      anon_sym_DQUOTE,
    ACTIONS(528), 1,
      anon_sym_for,
    ACTIONS(530), 1,
      anon_sym_LBRACE,
    STATE(108), 1,
      sym_block,
    STATE(226), 1,
      sym_string_literal,
    STATE(234), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4602] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(524), 1,
      aux_sym__pterm_token1,
    ACTIONS(526), 1,
      aux_sym__nterm_token1,
    ACTIONS(532), 3,
      anon_sym_DASH_GT,
      anon_sym_RBRACE,
      anon_sym_COLON,
    STATE(146), 3,
      sym__pterm,
      sym__nterm,
      aux_sym_chems_repeat1,
  [4622] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(260), 2,
      anon_sym_for,
      sym_identifier,
    ACTIONS(262), 6,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
      anon_sym_SQUOTE,
      anon_sym_DQUOTE,
      anon_sym_LBRACE,
  [4638] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(534), 1,
      sym_identifier,
    STATE(154), 1,
      sym_obj_key,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [4657] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(534), 1,
      sym_identifier,
    STATE(155), 1,
      sym_obj_key,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [4676] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(534), 1,
      sym_identifier,
    STATE(153), 1,
      sym_obj_key,
    STATE(141), 2,
      sym_special_keys,
      aux_sym_obj_key_repeat1,
    ACTIONS(19), 3,
      anon_sym_1d,
      anon_sym_2d,
      anon_sym_3d,
  [4695] = 7,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(236), 1,
      anon_sym_SQUOTE,
    ACTIONS(238), 1,
      anon_sym_DQUOTE,
    ACTIONS(512), 1,
      anon_sym_LBRACE,
    STATE(58), 1,
      sym_block,
    STATE(241), 1,
      sym_string_literal,
    STATE(77), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4718] = 7,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(236), 1,
      anon_sym_SQUOTE,
    ACTIONS(238), 1,
      anon_sym_DQUOTE,
    ACTIONS(536), 1,
      anon_sym_LBRACE,
    STATE(129), 1,
      sym_block,
    STATE(233), 1,
      sym_string_literal,
    STATE(77), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4741] = 7,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(236), 1,
      anon_sym_SQUOTE,
    ACTIONS(238), 1,
      anon_sym_DQUOTE,
    ACTIONS(530), 1,
      anon_sym_LBRACE,
    STATE(112), 1,
      sym_block,
    STATE(244), 1,
      sym_string_literal,
    STATE(77), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4764] = 6,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(236), 1,
      anon_sym_SQUOTE,
    ACTIONS(238), 1,
      anon_sym_DQUOTE,
    ACTIONS(538), 1,
      sym_identifier,
    STATE(96), 1,
      sym_string_literal,
    STATE(77), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4784] = 6,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(252), 1,
      anon_sym_SQUOTE,
    ACTIONS(254), 1,
      anon_sym_DQUOTE,
    ACTIONS(540), 1,
      sym_identifier,
    STATE(122), 1,
      sym_string_literal,
    STATE(134), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4804] = 6,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(236), 1,
      anon_sym_SQUOTE,
    ACTIONS(238), 1,
      anon_sym_DQUOTE,
    ACTIONS(542), 1,
      sym_identifier,
    STATE(238), 1,
      sym_string_literal,
    STATE(77), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4824] = 6,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(236), 1,
      anon_sym_SQUOTE,
    ACTIONS(238), 1,
      anon_sym_DQUOTE,
    ACTIONS(544), 1,
      sym_identifier,
    STATE(250), 1,
      sym_string_literal,
    STATE(77), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4844] = 6,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(252), 1,
      anon_sym_SQUOTE,
    ACTIONS(254), 1,
      anon_sym_DQUOTE,
    ACTIONS(546), 1,
      sym_identifier,
    STATE(135), 1,
      sym_string_literal,
    STATE(134), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4864] = 6,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(236), 1,
      anon_sym_SQUOTE,
    ACTIONS(238), 1,
      anon_sym_DQUOTE,
    ACTIONS(548), 1,
      sym_identifier,
    STATE(246), 1,
      sym_string_literal,
    STATE(77), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4884] = 6,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(236), 1,
      anon_sym_SQUOTE,
    ACTIONS(238), 1,
      anon_sym_DQUOTE,
    ACTIONS(550), 1,
      sym_identifier,
    STATE(248), 1,
      sym_string_literal,
    STATE(77), 2,
      sym__string_literal_single,
      sym__string_literal_double,
  [4904] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(552), 5,
      anon_sym_DASH_GT,
      anon_sym_RBRACE,
      aux_sym__pterm_token1,
      aux_sym__nterm_token1,
      anon_sym_COLON,
  [4915] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(554), 1,
      anon_sym_rate,
    ACTIONS(556), 1,
      anon_sym_janev_rate,
    STATE(40), 1,
      sym_rate,
    STATE(61), 2,
      sym_arrhenius,
      sym_janev,
  [4932] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(558), 1,
      anon_sym_rate,
    ACTIONS(560), 1,
      anon_sym_janev_rate,
    STATE(249), 1,
      sym_rate,
    STATE(294), 2,
      sym_arrhenius,
      sym_janev,
  [4949] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(558), 1,
      anon_sym_rate,
    ACTIONS(560), 1,
      anon_sym_janev_rate,
    STATE(247), 1,
      sym_rate,
    STATE(294), 2,
      sym_arrhenius,
      sym_janev,
  [4966] = 5,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(562), 1,
      anon_sym_rate,
    ACTIONS(564), 1,
      anon_sym_janev_rate,
    STATE(89), 1,
      sym_rate,
    STATE(92), 2,
      sym_arrhenius,
      sym_janev,
  [4983] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(566), 5,
      anon_sym_DASH_GT,
      anon_sym_RBRACE,
      aux_sym__pterm_token1,
      aux_sym__nterm_token1,
      anon_sym_COLON,
  [4994] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(570), 1,
      anon_sym_RPAREN,
    STATE(267), 1,
      sym__rawqty,
    ACTIONS(568), 2,
      sym_define_ref,
      sym_decimal,
  [5008] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(574), 1,
      anon_sym_RPAREN,
    STATE(251), 1,
      sym__rawqty,
    ACTIONS(572), 2,
      sym_define_ref,
      sym_decimal,
  [5022] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(578), 1,
      anon_sym_COLON,
    STATE(287), 1,
      sym__rawqty,
    ACTIONS(576), 2,
      sym_define_ref,
      sym_decimal,
  [5036] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(582), 1,
      anon_sym_RPAREN,
    STATE(257), 1,
      sym__rawqty,
    ACTIONS(580), 2,
      sym_define_ref,
      sym_decimal,
  [5050] = 4,
    ACTIONS(3), 1,
      sym_comment,
    STATE(232), 1,
      sym_sub_formula,
    STATE(308), 1,
      sym_chems,
    ACTIONS(584), 2,
      sym_define_key,
      sym_identifier,
  [5064] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(588), 1,
      anon_sym_COLON,
    STATE(320), 1,
      sym__rawqty,
    ACTIONS(586), 2,
      sym_define_ref,
      sym_decimal,
  [5078] = 4,
    ACTIONS(3), 1,
      sym_comment,
    STATE(240), 1,
      sym_sub_formula,
    STATE(308), 1,
      sym_chems,
    ACTIONS(584), 2,
      sym_define_key,
      sym_identifier,
  [5092] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(163), 1,
      sym__rawqty,
    ACTIONS(590), 3,
      sym_define_ref,
      sym_decimal,
      sym_identifier,
  [5104] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(594), 1,
      anon_sym_RPAREN,
    STATE(329), 1,
      sym__rawqty,
    ACTIONS(592), 2,
      sym_define_ref,
      sym_decimal,
  [5118] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(216), 1,
      sym__rawqty,
    ACTIONS(596), 2,
      sym_define_ref,
      sym_decimal,
  [5129] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(512), 1,
      anon_sym_LBRACE,
    ACTIONS(598), 1,
      anon_sym_for,
    STATE(51), 1,
      sym_block,
  [5142] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(186), 1,
      sym__rawqty,
    ACTIONS(600), 2,
      sym_define_ref,
      sym_decimal,
  [5153] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(222), 1,
      sym__rawqty,
    ACTIONS(602), 2,
      sym_define_ref,
      sym_decimal,
  [5164] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(334), 1,
      sym__rawqty,
    ACTIONS(604), 2,
      sym_define_ref,
      sym_decimal,
  [5175] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(181), 1,
      sym__rawqty,
    ACTIONS(606), 2,
      sym_define_ref,
      sym_decimal,
  [5186] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(192), 1,
      sym__rawqty,
    ACTIONS(608), 2,
      sym_define_ref,
      sym_decimal,
  [5197] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(183), 1,
      sym__rawqty,
    ACTIONS(610), 2,
      sym_define_ref,
      sym_decimal,
  [5208] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(194), 1,
      sym__rawqty,
    ACTIONS(612), 2,
      sym_define_ref,
      sym_decimal,
  [5219] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(325), 1,
      sym__rawqty,
    ACTIONS(614), 2,
      sym_define_ref,
      sym_decimal,
  [5230] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(185), 1,
      sym__rawqty,
    ACTIONS(616), 2,
      sym_define_ref,
      sym_decimal,
  [5241] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(330), 1,
      sym__rawqty,
    ACTIONS(618), 2,
      sym_define_ref,
      sym_decimal,
  [5252] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(91), 1,
      sym__rawqty,
    ACTIONS(620), 2,
      sym_define_ref,
      sym_decimal,
  [5263] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(81), 1,
      sym__rawqty,
    ACTIONS(622), 2,
      sym_define_ref,
      sym_decimal,
  [5274] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(196), 1,
      sym__rawqty,
    ACTIONS(624), 2,
      sym_define_ref,
      sym_decimal,
  [5285] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(188), 1,
      sym__rawqty,
    ACTIONS(626), 2,
      sym_define_ref,
      sym_decimal,
  [5296] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(198), 1,
      sym__rawqty,
    ACTIONS(628), 2,
      sym_define_ref,
      sym_decimal,
  [5307] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(167), 1,
      sym__rawqty,
    ACTIONS(630), 2,
      sym_define_ref,
      sym_decimal,
  [5318] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(199), 1,
      sym__rawqty,
    ACTIONS(632), 2,
      sym_define_ref,
      sym_decimal,
  [5329] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(111), 1,
      sym__rawqty,
    ACTIONS(634), 2,
      sym_define_ref,
      sym_decimal,
  [5340] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(201), 1,
      sym__rawqty,
    ACTIONS(636), 2,
      sym_define_ref,
      sym_decimal,
  [5351] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(204), 1,
      sym__rawqty,
    ACTIONS(638), 2,
      sym_define_ref,
      sym_decimal,
  [5362] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(193), 1,
      sym__rawqty,
    ACTIONS(640), 2,
      sym_define_ref,
      sym_decimal,
  [5373] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(205), 1,
      sym__rawqty,
    ACTIONS(642), 2,
      sym_define_ref,
      sym_decimal,
  [5384] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(326), 1,
      sym__rawqty,
    ACTIONS(644), 2,
      sym_define_ref,
      sym_decimal,
  [5395] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(309), 1,
      sym__rawqty,
    ACTIONS(646), 2,
      sym_define_ref,
      sym_decimal,
  [5406] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(206), 1,
      sym__rawqty,
    ACTIONS(648), 2,
      sym_define_ref,
      sym_decimal,
  [5417] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(178), 1,
      sym__rawqty,
    ACTIONS(650), 2,
      sym_define_ref,
      sym_decimal,
  [5428] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(212), 1,
      sym__rawqty,
    ACTIONS(652), 2,
      sym_define_ref,
      sym_decimal,
  [5439] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(197), 1,
      sym__rawqty,
    ACTIONS(654), 2,
      sym_define_ref,
      sym_decimal,
  [5450] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(200), 1,
      sym__rawqty,
    ACTIONS(656), 2,
      sym_define_ref,
      sym_decimal,
  [5461] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(202), 1,
      sym__rawqty,
    ACTIONS(658), 2,
      sym_define_ref,
      sym_decimal,
  [5472] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(660), 1,
      anon_sym_RBRACE,
    ACTIONS(662), 1,
      anon_sym_COLON,
    STATE(210), 1,
      aux_sym_full_formula_repeat1,
  [5485] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(217), 1,
      sym__rawqty,
    ACTIONS(665), 2,
      sym_define_ref,
      sym_decimal,
  [5496] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(207), 1,
      sym__rawqty,
    ACTIONS(667), 2,
      sym_define_ref,
      sym_decimal,
  [5507] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(164), 1,
      sym__rawqty,
    ACTIONS(669), 2,
      sym_define_ref,
      sym_decimal,
  [5518] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(228), 1,
      sym__rawqty,
    ACTIONS(671), 2,
      sym_define_ref,
      sym_decimal,
  [5529] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(673), 1,
      anon_sym_coulomb,
    ACTIONS(675), 1,
      anon_sym_cross,
    ACTIONS(677), 1,
      anon_sym_metallic,
  [5542] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(218), 1,
      sym__rawqty,
    ACTIONS(679), 2,
      sym_define_ref,
      sym_decimal,
  [5553] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(225), 1,
      sym__rawqty,
    ACTIONS(681), 2,
      sym_define_ref,
      sym_decimal,
  [5564] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(67), 1,
      sym__rawqty,
    ACTIONS(683), 2,
      sym_define_ref,
      sym_decimal,
  [5575] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(208), 1,
      sym__rawqty,
    ACTIONS(685), 2,
      sym_define_ref,
      sym_decimal,
  [5586] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(209), 1,
      sym__rawqty,
    ACTIONS(687), 2,
      sym_define_ref,
      sym_decimal,
  [5597] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(689), 1,
      anon_sym_RBRACE,
    ACTIONS(691), 1,
      anon_sym_COLON,
    STATE(210), 1,
      aux_sym_full_formula_repeat1,
  [5610] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(289), 1,
      sym__rawqty,
    ACTIONS(693), 2,
      sym_define_ref,
      sym_decimal,
  [5621] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(242), 1,
      sym_chems,
    ACTIONS(584), 2,
      sym_define_key,
      sym_identifier,
  [5632] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(168), 1,
      sym__rawqty,
    ACTIONS(695), 2,
      sym_define_ref,
      sym_decimal,
  [5643] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(64), 1,
      sym__rawqty,
    ACTIONS(697), 2,
      sym_define_ref,
      sym_decimal,
  [5654] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(530), 1,
      anon_sym_LBRACE,
    ACTIONS(699), 1,
      anon_sym_for,
    STATE(97), 1,
      sym_block,
  [5667] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(39), 1,
      sym__rawqty,
    ACTIONS(701), 2,
      sym_define_ref,
      sym_decimal,
  [5678] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(191), 1,
      sym__rawqty,
    ACTIONS(703), 2,
      sym_define_ref,
      sym_decimal,
  [5689] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(705), 1,
      anon_sym_coulomb,
    ACTIONS(707), 1,
      anon_sym_cross,
    ACTIONS(709), 1,
      anon_sym_metallic,
  [5702] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(100), 1,
      sym__rawqty,
    ACTIONS(711), 2,
      sym_define_ref,
      sym_decimal,
  [5713] = 3,
    ACTIONS(3), 1,
      sym_comment,
    STATE(46), 1,
      sym__rawqty,
    ACTIONS(713), 2,
      sym_define_ref,
      sym_decimal,
  [5724] = 4,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(691), 1,
      anon_sym_COLON,
    ACTIONS(715), 1,
      anon_sym_RBRACE,
    STATE(221), 1,
      aux_sym_full_formula_repeat1,
  [5737] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(536), 1,
      anon_sym_LBRACE,
    STATE(123), 1,
      sym_block,
  [5747] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(322), 2,
      anon_sym_for,
      anon_sym_LBRACE,
  [5755] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(314), 2,
      anon_sym_for,
      anon_sym_LBRACE,
  [5763] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(717), 1,
      anon_sym_LBRACE,
    STATE(166), 1,
      sym_full_formula,
  [5773] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(350), 2,
      anon_sym_for,
      anon_sym_LBRACE,
  [5781] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(530), 1,
      anon_sym_LBRACE,
    STATE(90), 1,
      sym_block,
  [5791] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(719), 2,
      anon_sym_rate,
      anon_sym_janev_rate,
  [5799] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(660), 2,
      anon_sym_RBRACE,
      anon_sym_COLON,
  [5807] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(512), 1,
      anon_sym_LBRACE,
    STATE(66), 1,
      sym_block,
  [5817] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(721), 2,
      anon_sym_RBRACE,
      anon_sym_COLON,
  [5825] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(717), 1,
      anon_sym_LBRACE,
    STATE(165), 1,
      sym_full_formula,
  [5835] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(530), 1,
      anon_sym_LBRACE,
    STATE(87), 1,
      sym_block,
  [5845] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(723), 2,
      anon_sym_rate,
      anon_sym_janev_rate,
  [5853] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(512), 1,
      anon_sym_LBRACE,
    STATE(69), 1,
      sym_block,
  [5863] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(725), 1,
      sym_identifier,
    STATE(68), 1,
      sym_range,
  [5873] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(512), 1,
      anon_sym_LBRACE,
    STATE(63), 1,
      sym_block,
  [5883] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(727), 1,
      sym_identifier,
    STATE(102), 1,
      sym_range,
  [5893] = 3,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(530), 1,
      anon_sym_LBRACE,
    STATE(105), 1,
      sym_block,
  [5903] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(594), 1,
      anon_sym_RPAREN,
  [5910] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(729), 1,
      ts_builtin_sym_end,
  [5917] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(731), 1,
      anon_sym_POUNDendif,
  [5924] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(733), 1,
      anon_sym_POUNDendif,
  [5931] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(735), 1,
      sym_identifier,
  [5938] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(737), 1,
      anon_sym_EQ,
  [5945] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(570), 1,
      anon_sym_RPAREN,
  [5952] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(739), 1,
      sym_identifier,
  [5959] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(741), 1,
      sym_identifier,
  [5966] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(743), 1,
      sym_identifier,
  [5973] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(745), 1,
      anon_sym_SQUOTE,
  [5980] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(747), 1,
      anon_sym_DQUOTE,
  [5987] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(749), 1,
      sym_identifier,
  [5994] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(134), 1,
      anon_sym_POUNDendif,
  [6001] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(751), 1,
      anon_sym_DQUOTE,
  [6008] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(753), 1,
      anon_sym_LPAREN,
  [6015] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(755), 1,
      anon_sym_RPAREN,
  [6022] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(757), 1,
      anon_sym_SQUOTE,
  [6029] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(759), 1,
      anon_sym_DQUOTE,
  [6036] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(761), 1,
      anon_sym_SQUOTE,
  [6043] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(763), 1,
      anon_sym_DQUOTE,
  [6050] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(765), 1,
      sym_define_key,
  [6057] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(767), 1,
      sym_define_key,
  [6064] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(769), 1,
      anon_sym_SQUOTE,
  [6071] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(771), 1,
      anon_sym_EQ,
  [6078] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(773), 1,
      anon_sym_EQ,
  [6085] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(775), 1,
      sym_identifier,
  [6092] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(777), 1,
      sym_identifier,
  [6099] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(779), 1,
      anon_sym_EQ,
  [6106] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(781), 1,
      anon_sym_EQ,
  [6113] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(783), 1,
      sym_identifier,
  [6120] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(785), 1,
      sym_identifier,
  [6127] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(787), 1,
      anon_sym_EQ,
  [6134] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(789), 1,
      anon_sym_EQ,
  [6141] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(791), 1,
      anon_sym_EQ,
  [6148] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(793), 1,
      anon_sym_SQUOTE,
  [6155] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(795), 1,
      anon_sym_COLON,
  [6162] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(797), 1,
      anon_sym_EQ,
  [6169] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(412), 1,
      sym_identifier,
  [6176] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(799), 1,
      anon_sym_LT_DASH_GT,
  [6183] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(801), 1,
      anon_sym_EQ,
  [6190] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(803), 1,
      sym_identifier,
  [6197] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(805), 1,
      sym_identifier,
  [6204] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(390), 1,
      sym_identifier,
  [6211] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(807), 1,
      sym_define_key,
  [6218] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(809), 1,
      anon_sym_EQ,
  [6225] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(811), 1,
      anon_sym_level,
  [6232] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(813), 1,
      anon_sym_DASH_GT,
  [6239] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(815), 1,
      sym_identifier,
  [6246] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(817), 1,
      sym_identifier,
  [6253] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(819), 1,
      sym_identifier,
  [6260] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(128), 1,
      anon_sym_POUNDendif,
  [6267] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(821), 1,
      anon_sym_LT_DASH_GT,
  [6274] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(823), 1,
      anon_sym_EQ,
  [6281] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(825), 1,
      anon_sym_LPAREN,
  [6288] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(827), 1,
      anon_sym_section,
  [6295] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(829), 1,
      anon_sym_EQ,
  [6302] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(831), 1,
      anon_sym_DASH_GT,
  [6309] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(833), 1,
      anon_sym_ref_density,
  [6316] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(835), 1,
      anon_sym_EQ,
  [6323] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(837), 1,
      anon_sym_section,
  [6330] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(839), 1,
      anon_sym_EQ,
  [6337] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(841), 1,
      sym_identifier,
  [6344] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(843), 1,
      anon_sym_level,
  [6351] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(845), 1,
      anon_sym_ks,
  [6358] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(847), 1,
      anon_sym_EQ,
  [6365] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(849), 1,
      anon_sym_DQUOTE,
  [6372] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(851), 1,
      anon_sym_EQ,
  [6379] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(853), 1,
      sym_identifier,
  [6386] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(855), 1,
      anon_sym_COLON,
  [6393] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(857), 1,
      anon_sym_EQ,
  [6400] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(859), 1,
      anon_sym_EQ,
  [6407] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(861), 1,
      anon_sym_DASH_GT,
  [6414] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(863), 1,
      anon_sym_EQ,
  [6421] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(865), 1,
      anon_sym_fermi_energy_ev,
  [6428] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(402), 1,
      sym_identifier,
  [6435] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(867), 1,
      sym_identifier,
  [6442] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(869), 1,
      sym_define_key,
  [6449] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(871), 1,
      anon_sym_RPAREN,
  [6456] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(873), 1,
      anon_sym_fermi_energy_ev,
  [6463] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(875), 1,
      anon_sym_EQ,
  [6470] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(877), 1,
      sym_identifier,
  [6477] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(879), 1,
      anon_sym_EQ,
  [6484] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(881), 1,
      anon_sym_ref_density,
  [6491] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(883), 1,
      anon_sym_EQ,
  [6498] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(885), 1,
      anon_sym_ks,
  [6505] = 2,
    ACTIONS(3), 1,
      sym_comment,
    ACTIONS(887), 1,
      anon_sym_EQ,
};

static const uint32_t ts_small_parse_table_map[] = {
  [SMALL_STATE(2)] = 0,
  [SMALL_STATE(3)] = 69,
  [SMALL_STATE(4)] = 132,
  [SMALL_STATE(5)] = 195,
  [SMALL_STATE(6)] = 256,
  [SMALL_STATE(7)] = 336,
  [SMALL_STATE(8)] = 416,
  [SMALL_STATE(9)] = 493,
  [SMALL_STATE(10)] = 574,
  [SMALL_STATE(11)] = 651,
  [SMALL_STATE(12)] = 728,
  [SMALL_STATE(13)] = 788,
  [SMALL_STATE(14)] = 848,
  [SMALL_STATE(15)] = 908,
  [SMALL_STATE(16)] = 964,
  [SMALL_STATE(17)] = 1024,
  [SMALL_STATE(18)] = 1078,
  [SMALL_STATE(19)] = 1132,
  [SMALL_STATE(20)] = 1186,
  [SMALL_STATE(21)] = 1236,
  [SMALL_STATE(22)] = 1286,
  [SMALL_STATE(23)] = 1336,
  [SMALL_STATE(24)] = 1386,
  [SMALL_STATE(25)] = 1436,
  [SMALL_STATE(26)] = 1486,
  [SMALL_STATE(27)] = 1536,
  [SMALL_STATE(28)] = 1572,
  [SMALL_STATE(29)] = 1604,
  [SMALL_STATE(30)] = 1649,
  [SMALL_STATE(31)] = 1694,
  [SMALL_STATE(32)] = 1739,
  [SMALL_STATE(33)] = 1767,
  [SMALL_STATE(34)] = 1797,
  [SMALL_STATE(35)] = 1827,
  [SMALL_STATE(36)] = 1869,
  [SMALL_STATE(37)] = 1896,
  [SMALL_STATE(38)] = 1923,
  [SMALL_STATE(39)] = 1950,
  [SMALL_STATE(40)] = 1977,
  [SMALL_STATE(41)] = 2004,
  [SMALL_STATE(42)] = 2031,
  [SMALL_STATE(43)] = 2058,
  [SMALL_STATE(44)] = 2085,
  [SMALL_STATE(45)] = 2126,
  [SMALL_STATE(46)] = 2167,
  [SMALL_STATE(47)] = 2194,
  [SMALL_STATE(48)] = 2235,
  [SMALL_STATE(49)] = 2262,
  [SMALL_STATE(50)] = 2289,
  [SMALL_STATE(51)] = 2316,
  [SMALL_STATE(52)] = 2343,
  [SMALL_STATE(53)] = 2370,
  [SMALL_STATE(54)] = 2397,
  [SMALL_STATE(55)] = 2424,
  [SMALL_STATE(56)] = 2451,
  [SMALL_STATE(57)] = 2478,
  [SMALL_STATE(58)] = 2505,
  [SMALL_STATE(59)] = 2532,
  [SMALL_STATE(60)] = 2573,
  [SMALL_STATE(61)] = 2614,
  [SMALL_STATE(62)] = 2641,
  [SMALL_STATE(63)] = 2668,
  [SMALL_STATE(64)] = 2695,
  [SMALL_STATE(65)] = 2722,
  [SMALL_STATE(66)] = 2763,
  [SMALL_STATE(67)] = 2790,
  [SMALL_STATE(68)] = 2817,
  [SMALL_STATE(69)] = 2844,
  [SMALL_STATE(70)] = 2871,
  [SMALL_STATE(71)] = 2909,
  [SMALL_STATE(72)] = 2947,
  [SMALL_STATE(73)] = 2985,
  [SMALL_STATE(74)] = 3023,
  [SMALL_STATE(75)] = 3061,
  [SMALL_STATE(76)] = 3099,
  [SMALL_STATE(77)] = 3125,
  [SMALL_STATE(78)] = 3147,
  [SMALL_STATE(79)] = 3169,
  [SMALL_STATE(80)] = 3191,
  [SMALL_STATE(81)] = 3212,
  [SMALL_STATE(82)] = 3233,
  [SMALL_STATE(83)] = 3254,
  [SMALL_STATE(84)] = 3275,
  [SMALL_STATE(85)] = 3296,
  [SMALL_STATE(86)] = 3317,
  [SMALL_STATE(87)] = 3338,
  [SMALL_STATE(88)] = 3359,
  [SMALL_STATE(89)] = 3380,
  [SMALL_STATE(90)] = 3401,
  [SMALL_STATE(91)] = 3422,
  [SMALL_STATE(92)] = 3443,
  [SMALL_STATE(93)] = 3464,
  [SMALL_STATE(94)] = 3485,
  [SMALL_STATE(95)] = 3506,
  [SMALL_STATE(96)] = 3527,
  [SMALL_STATE(97)] = 3548,
  [SMALL_STATE(98)] = 3569,
  [SMALL_STATE(99)] = 3590,
  [SMALL_STATE(100)] = 3611,
  [SMALL_STATE(101)] = 3632,
  [SMALL_STATE(102)] = 3653,
  [SMALL_STATE(103)] = 3674,
  [SMALL_STATE(104)] = 3695,
  [SMALL_STATE(105)] = 3716,
  [SMALL_STATE(106)] = 3737,
  [SMALL_STATE(107)] = 3758,
  [SMALL_STATE(108)] = 3779,
  [SMALL_STATE(109)] = 3800,
  [SMALL_STATE(110)] = 3821,
  [SMALL_STATE(111)] = 3842,
  [SMALL_STATE(112)] = 3863,
  [SMALL_STATE(113)] = 3884,
  [SMALL_STATE(114)] = 3905,
  [SMALL_STATE(115)] = 3926,
  [SMALL_STATE(116)] = 3947,
  [SMALL_STATE(117)] = 3968,
  [SMALL_STATE(118)] = 3987,
  [SMALL_STATE(119)] = 4006,
  [SMALL_STATE(120)] = 4025,
  [SMALL_STATE(121)] = 4044,
  [SMALL_STATE(122)] = 4063,
  [SMALL_STATE(123)] = 4082,
  [SMALL_STATE(124)] = 4101,
  [SMALL_STATE(125)] = 4120,
  [SMALL_STATE(126)] = 4139,
  [SMALL_STATE(127)] = 4158,
  [SMALL_STATE(128)] = 4177,
  [SMALL_STATE(129)] = 4196,
  [SMALL_STATE(130)] = 4215,
  [SMALL_STATE(131)] = 4234,
  [SMALL_STATE(132)] = 4253,
  [SMALL_STATE(133)] = 4272,
  [SMALL_STATE(134)] = 4291,
  [SMALL_STATE(135)] = 4310,
  [SMALL_STATE(136)] = 4329,
  [SMALL_STATE(137)] = 4348,
  [SMALL_STATE(138)] = 4372,
  [SMALL_STATE(139)] = 4400,
  [SMALL_STATE(140)] = 4422,
  [SMALL_STATE(141)] = 4446,
  [SMALL_STATE(142)] = 4468,
  [SMALL_STATE(143)] = 4496,
  [SMALL_STATE(144)] = 4522,
  [SMALL_STATE(145)] = 4536,
  [SMALL_STATE(146)] = 4556,
  [SMALL_STATE(147)] = 4576,
  [SMALL_STATE(148)] = 4602,
  [SMALL_STATE(149)] = 4622,
  [SMALL_STATE(150)] = 4638,
  [SMALL_STATE(151)] = 4657,
  [SMALL_STATE(152)] = 4676,
  [SMALL_STATE(153)] = 4695,
  [SMALL_STATE(154)] = 4718,
  [SMALL_STATE(155)] = 4741,
  [SMALL_STATE(156)] = 4764,
  [SMALL_STATE(157)] = 4784,
  [SMALL_STATE(158)] = 4804,
  [SMALL_STATE(159)] = 4824,
  [SMALL_STATE(160)] = 4844,
  [SMALL_STATE(161)] = 4864,
  [SMALL_STATE(162)] = 4884,
  [SMALL_STATE(163)] = 4904,
  [SMALL_STATE(164)] = 4915,
  [SMALL_STATE(165)] = 4932,
  [SMALL_STATE(166)] = 4949,
  [SMALL_STATE(167)] = 4966,
  [SMALL_STATE(168)] = 4983,
  [SMALL_STATE(169)] = 4994,
  [SMALL_STATE(170)] = 5008,
  [SMALL_STATE(171)] = 5022,
  [SMALL_STATE(172)] = 5036,
  [SMALL_STATE(173)] = 5050,
  [SMALL_STATE(174)] = 5064,
  [SMALL_STATE(175)] = 5078,
  [SMALL_STATE(176)] = 5092,
  [SMALL_STATE(177)] = 5104,
  [SMALL_STATE(178)] = 5118,
  [SMALL_STATE(179)] = 5129,
  [SMALL_STATE(180)] = 5142,
  [SMALL_STATE(181)] = 5153,
  [SMALL_STATE(182)] = 5164,
  [SMALL_STATE(183)] = 5175,
  [SMALL_STATE(184)] = 5186,
  [SMALL_STATE(185)] = 5197,
  [SMALL_STATE(186)] = 5208,
  [SMALL_STATE(187)] = 5219,
  [SMALL_STATE(188)] = 5230,
  [SMALL_STATE(189)] = 5241,
  [SMALL_STATE(190)] = 5252,
  [SMALL_STATE(191)] = 5263,
  [SMALL_STATE(192)] = 5274,
  [SMALL_STATE(193)] = 5285,
  [SMALL_STATE(194)] = 5296,
  [SMALL_STATE(195)] = 5307,
  [SMALL_STATE(196)] = 5318,
  [SMALL_STATE(197)] = 5329,
  [SMALL_STATE(198)] = 5340,
  [SMALL_STATE(199)] = 5351,
  [SMALL_STATE(200)] = 5362,
  [SMALL_STATE(201)] = 5373,
  [SMALL_STATE(202)] = 5384,
  [SMALL_STATE(203)] = 5395,
  [SMALL_STATE(204)] = 5406,
  [SMALL_STATE(205)] = 5417,
  [SMALL_STATE(206)] = 5428,
  [SMALL_STATE(207)] = 5439,
  [SMALL_STATE(208)] = 5450,
  [SMALL_STATE(209)] = 5461,
  [SMALL_STATE(210)] = 5472,
  [SMALL_STATE(211)] = 5485,
  [SMALL_STATE(212)] = 5496,
  [SMALL_STATE(213)] = 5507,
  [SMALL_STATE(214)] = 5518,
  [SMALL_STATE(215)] = 5529,
  [SMALL_STATE(216)] = 5542,
  [SMALL_STATE(217)] = 5553,
  [SMALL_STATE(218)] = 5564,
  [SMALL_STATE(219)] = 5575,
  [SMALL_STATE(220)] = 5586,
  [SMALL_STATE(221)] = 5597,
  [SMALL_STATE(222)] = 5610,
  [SMALL_STATE(223)] = 5621,
  [SMALL_STATE(224)] = 5632,
  [SMALL_STATE(225)] = 5643,
  [SMALL_STATE(226)] = 5654,
  [SMALL_STATE(227)] = 5667,
  [SMALL_STATE(228)] = 5678,
  [SMALL_STATE(229)] = 5689,
  [SMALL_STATE(230)] = 5702,
  [SMALL_STATE(231)] = 5713,
  [SMALL_STATE(232)] = 5724,
  [SMALL_STATE(233)] = 5737,
  [SMALL_STATE(234)] = 5747,
  [SMALL_STATE(235)] = 5755,
  [SMALL_STATE(236)] = 5763,
  [SMALL_STATE(237)] = 5773,
  [SMALL_STATE(238)] = 5781,
  [SMALL_STATE(239)] = 5791,
  [SMALL_STATE(240)] = 5799,
  [SMALL_STATE(241)] = 5807,
  [SMALL_STATE(242)] = 5817,
  [SMALL_STATE(243)] = 5825,
  [SMALL_STATE(244)] = 5835,
  [SMALL_STATE(245)] = 5845,
  [SMALL_STATE(246)] = 5853,
  [SMALL_STATE(247)] = 5863,
  [SMALL_STATE(248)] = 5873,
  [SMALL_STATE(249)] = 5883,
  [SMALL_STATE(250)] = 5893,
  [SMALL_STATE(251)] = 5903,
  [SMALL_STATE(252)] = 5910,
  [SMALL_STATE(253)] = 5917,
  [SMALL_STATE(254)] = 5924,
  [SMALL_STATE(255)] = 5931,
  [SMALL_STATE(256)] = 5938,
  [SMALL_STATE(257)] = 5945,
  [SMALL_STATE(258)] = 5952,
  [SMALL_STATE(259)] = 5959,
  [SMALL_STATE(260)] = 5966,
  [SMALL_STATE(261)] = 5973,
  [SMALL_STATE(262)] = 5980,
  [SMALL_STATE(263)] = 5987,
  [SMALL_STATE(264)] = 5994,
  [SMALL_STATE(265)] = 6001,
  [SMALL_STATE(266)] = 6008,
  [SMALL_STATE(267)] = 6015,
  [SMALL_STATE(268)] = 6022,
  [SMALL_STATE(269)] = 6029,
  [SMALL_STATE(270)] = 6036,
  [SMALL_STATE(271)] = 6043,
  [SMALL_STATE(272)] = 6050,
  [SMALL_STATE(273)] = 6057,
  [SMALL_STATE(274)] = 6064,
  [SMALL_STATE(275)] = 6071,
  [SMALL_STATE(276)] = 6078,
  [SMALL_STATE(277)] = 6085,
  [SMALL_STATE(278)] = 6092,
  [SMALL_STATE(279)] = 6099,
  [SMALL_STATE(280)] = 6106,
  [SMALL_STATE(281)] = 6113,
  [SMALL_STATE(282)] = 6120,
  [SMALL_STATE(283)] = 6127,
  [SMALL_STATE(284)] = 6134,
  [SMALL_STATE(285)] = 6141,
  [SMALL_STATE(286)] = 6148,
  [SMALL_STATE(287)] = 6155,
  [SMALL_STATE(288)] = 6162,
  [SMALL_STATE(289)] = 6169,
  [SMALL_STATE(290)] = 6176,
  [SMALL_STATE(291)] = 6183,
  [SMALL_STATE(292)] = 6190,
  [SMALL_STATE(293)] = 6197,
  [SMALL_STATE(294)] = 6204,
  [SMALL_STATE(295)] = 6211,
  [SMALL_STATE(296)] = 6218,
  [SMALL_STATE(297)] = 6225,
  [SMALL_STATE(298)] = 6232,
  [SMALL_STATE(299)] = 6239,
  [SMALL_STATE(300)] = 6246,
  [SMALL_STATE(301)] = 6253,
  [SMALL_STATE(302)] = 6260,
  [SMALL_STATE(303)] = 6267,
  [SMALL_STATE(304)] = 6274,
  [SMALL_STATE(305)] = 6281,
  [SMALL_STATE(306)] = 6288,
  [SMALL_STATE(307)] = 6295,
  [SMALL_STATE(308)] = 6302,
  [SMALL_STATE(309)] = 6309,
  [SMALL_STATE(310)] = 6316,
  [SMALL_STATE(311)] = 6323,
  [SMALL_STATE(312)] = 6330,
  [SMALL_STATE(313)] = 6337,
  [SMALL_STATE(314)] = 6344,
  [SMALL_STATE(315)] = 6351,
  [SMALL_STATE(316)] = 6358,
  [SMALL_STATE(317)] = 6365,
  [SMALL_STATE(318)] = 6372,
  [SMALL_STATE(319)] = 6379,
  [SMALL_STATE(320)] = 6386,
  [SMALL_STATE(321)] = 6393,
  [SMALL_STATE(322)] = 6400,
  [SMALL_STATE(323)] = 6407,
  [SMALL_STATE(324)] = 6414,
  [SMALL_STATE(325)] = 6421,
  [SMALL_STATE(326)] = 6428,
  [SMALL_STATE(327)] = 6435,
  [SMALL_STATE(328)] = 6442,
  [SMALL_STATE(329)] = 6449,
  [SMALL_STATE(330)] = 6456,
  [SMALL_STATE(331)] = 6463,
  [SMALL_STATE(332)] = 6470,
  [SMALL_STATE(333)] = 6477,
  [SMALL_STATE(334)] = 6484,
  [SMALL_STATE(335)] = 6491,
  [SMALL_STATE(336)] = 6498,
  [SMALL_STATE(337)] = 6505,
};

static const TSParseActionEntry ts_parse_actions[] = {
  [0] = {.entry = {.count = 0, .reusable = false}},
  [1] = {.entry = {.count = 1, .reusable = false}}, RECOVER(),
  [3] = {.entry = {.count = 1, .reusable = true}}, SHIFT_EXTRA(),
  [5] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_input_file, 0),
  [7] = {.entry = {.count = 1, .reusable = true}}, SHIFT(156),
  [9] = {.entry = {.count = 1, .reusable = true}}, SHIFT(328),
  [11] = {.entry = {.count = 1, .reusable = true}}, SHIFT(295),
  [13] = {.entry = {.count = 1, .reusable = false}}, SHIFT(138),
  [15] = {.entry = {.count = 1, .reusable = false}}, SHIFT(151),
  [17] = {.entry = {.count = 1, .reusable = false}}, SHIFT(141),
  [19] = {.entry = {.count = 1, .reusable = true}}, SHIFT(144),
  [21] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym__value, 1),
  [23] = {.entry = {.count = 1, .reusable = false}}, SHIFT(53),
  [25] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym__value, 1),
  [27] = {.entry = {.count = 1, .reusable = false}}, SHIFT_EXTRA(),
  [29] = {.entry = {.count = 1, .reusable = true}}, SHIFT(110),
  [31] = {.entry = {.count = 1, .reusable = true}}, SHIFT(99),
  [33] = {.entry = {.count = 1, .reusable = true}}, SHIFT(132),
  [35] = {.entry = {.count = 1, .reusable = false}}, SHIFT(142),
  [37] = {.entry = {.count = 1, .reusable = false}}, SHIFT(152),
  [39] = {.entry = {.count = 1, .reusable = false}}, SHIFT(10),
  [41] = {.entry = {.count = 1, .reusable = false}}, SHIFT(2),
  [43] = {.entry = {.count = 1, .reusable = false}}, SHIFT(57),
  [45] = {.entry = {.count = 1, .reusable = false}}, SHIFT(33),
  [47] = {.entry = {.count = 1, .reusable = false}}, SHIFT(32),
  [49] = {.entry = {.count = 1, .reusable = false}}, SHIFT(277),
  [51] = {.entry = {.count = 1, .reusable = false}}, SHIFT(278),
  [53] = {.entry = {.count = 1, .reusable = false}}, SHIFT(9),
  [55] = {.entry = {.count = 1, .reusable = false}}, SHIFT(74),
  [57] = {.entry = {.count = 1, .reusable = true}}, SHIFT(133),
  [59] = {.entry = {.count = 1, .reusable = false}}, SHIFT(8),
  [61] = {.entry = {.count = 1, .reusable = true}}, SHIFT(80),
  [63] = {.entry = {.count = 1, .reusable = false}}, SHIFT(11),
  [65] = {.entry = {.count = 1, .reusable = true}}, SHIFT(115),
  [67] = {.entry = {.count = 1, .reusable = true}}, SHIFT(160),
  [69] = {.entry = {.count = 1, .reusable = true}}, SHIFT(272),
  [71] = {.entry = {.count = 1, .reusable = true}}, SHIFT(273),
  [73] = {.entry = {.count = 1, .reusable = false}}, SHIFT(150),
  [75] = {.entry = {.count = 1, .reusable = false}}, SHIFT(157),
  [77] = {.entry = {.count = 1, .reusable = true}}, SHIFT(59),
  [79] = {.entry = {.count = 1, .reusable = false}}, SHIFT(3),
  [81] = {.entry = {.count = 1, .reusable = false}}, SHIFT(114),
  [83] = {.entry = {.count = 1, .reusable = false}}, SHIFT(76),
  [85] = {.entry = {.count = 1, .reusable = true}}, SHIFT(259),
  [87] = {.entry = {.count = 1, .reusable = true}}, SHIFT(299),
  [89] = {.entry = {.count = 1, .reusable = true}}, SHIFT(37),
  [91] = {.entry = {.count = 1, .reusable = true}}, SHIFT(125),
  [93] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(142),
  [96] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(152),
  [99] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(11),
  [102] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(2),
  [105] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(57),
  [108] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(33),
  [111] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(32),
  [114] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(277),
  [117] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(278),
  [120] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(9),
  [123] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 2), SHIFT_REPEAT(74),
  [126] = {.entry = {.count = 1, .reusable = true}}, REDUCE(aux_sym__define_value_repeat1, 2),
  [128] = {.entry = {.count = 1, .reusable = true}}, SHIFT(124),
  [130] = {.entry = {.count = 1, .reusable = true}}, SHIFT(19),
  [132] = {.entry = {.count = 1, .reusable = true}}, SHIFT(131),
  [134] = {.entry = {.count = 1, .reusable = true}}, SHIFT(101),
  [136] = {.entry = {.count = 1, .reusable = true}}, REDUCE(aux_sym_input_file_repeat1, 2),
  [138] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_input_file_repeat1, 2), SHIFT_REPEAT(156),
  [141] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_input_file_repeat1, 2), SHIFT_REPEAT(328),
  [144] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_input_file_repeat1, 2), SHIFT_REPEAT(295),
  [147] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_input_file_repeat1, 2), SHIFT_REPEAT(138),
  [150] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_input_file_repeat1, 2), SHIFT_REPEAT(151),
  [153] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_input_file_repeat1, 2), SHIFT_REPEAT(141),
  [156] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_input_file_repeat1, 2), SHIFT_REPEAT(144),
  [159] = {.entry = {.count = 1, .reusable = true}}, SHIFT(107),
  [161] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_else_block, 2),
  [163] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_input_file, 1),
  [165] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_else_block, 1),
  [167] = {.entry = {.count = 1, .reusable = true}}, SHIFT(52),
  [169] = {.entry = {.count = 1, .reusable = true}}, SHIFT(83),
  [171] = {.entry = {.count = 1, .reusable = true}}, SHIFT(121),
  [173] = {.entry = {.count = 1, .reusable = true}}, SHIFT(82),
  [175] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_block_repeat1, 2), SHIFT_REPEAT(160),
  [178] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_block_repeat1, 2), SHIFT_REPEAT(272),
  [181] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_block_repeat1, 2), SHIFT_REPEAT(273),
  [184] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_block_repeat1, 2), SHIFT_REPEAT(150),
  [187] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_block_repeat1, 2), SHIFT_REPEAT(157),
  [190] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_block_repeat1, 2), SHIFT_REPEAT(141),
  [193] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_block_repeat1, 2), SHIFT_REPEAT(144),
  [196] = {.entry = {.count = 1, .reusable = true}}, REDUCE(aux_sym_block_repeat1, 2),
  [198] = {.entry = {.count = 1, .reusable = true}}, SHIFT(126),
  [200] = {.entry = {.count = 1, .reusable = false}}, REDUCE(aux_sym_obj_key_repeat1, 2),
  [202] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_obj_key_repeat1, 2), SHIFT_REPEAT(27),
  [205] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_obj_key_repeat1, 2), SHIFT_REPEAT(32),
  [208] = {.entry = {.count = 1, .reusable = true}}, REDUCE(aux_sym_obj_key_repeat1, 2),
  [210] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_obj_key, 1),
  [212] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_obj_key, 1),
  [214] = {.entry = {.count = 1, .reusable = true}}, SHIFT(56),
  [216] = {.entry = {.count = 1, .reusable = true}}, SHIFT(2),
  [218] = {.entry = {.count = 1, .reusable = false}}, SHIFT(56),
  [220] = {.entry = {.count = 1, .reusable = true}}, SHIFT(277),
  [222] = {.entry = {.count = 1, .reusable = true}}, SHIFT(278),
  [224] = {.entry = {.count = 1, .reusable = true}}, SHIFT(70),
  [226] = {.entry = {.count = 1, .reusable = true}}, SHIFT(74),
  [228] = {.entry = {.count = 1, .reusable = true}}, SHIFT(116),
  [230] = {.entry = {.count = 1, .reusable = true}}, SHIFT(4),
  [232] = {.entry = {.count = 1, .reusable = false}}, SHIFT(103),
  [234] = {.entry = {.count = 1, .reusable = false}}, SHIFT(116),
  [236] = {.entry = {.count = 1, .reusable = true}}, SHIFT(282),
  [238] = {.entry = {.count = 1, .reusable = true}}, SHIFT(281),
  [240] = {.entry = {.count = 1, .reusable = true}}, SHIFT(72),
  [242] = {.entry = {.count = 1, .reusable = true}}, SHIFT(71),
  [244] = {.entry = {.count = 1, .reusable = true}}, SHIFT(120),
  [246] = {.entry = {.count = 1, .reusable = true}}, SHIFT(5),
  [248] = {.entry = {.count = 1, .reusable = false}}, SHIFT(119),
  [250] = {.entry = {.count = 1, .reusable = false}}, SHIFT(120),
  [252] = {.entry = {.count = 1, .reusable = true}}, SHIFT(300),
  [254] = {.entry = {.count = 1, .reusable = true}}, SHIFT(301),
  [256] = {.entry = {.count = 1, .reusable = true}}, SHIFT(75),
  [258] = {.entry = {.count = 1, .reusable = true}}, SHIFT(73),
  [260] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_special_keys, 1),
  [262] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_special_keys, 1),
  [264] = {.entry = {.count = 1, .reusable = false}}, REDUCE(aux_sym_obj_key_repeat1, 1),
  [266] = {.entry = {.count = 1, .reusable = false}}, SHIFT(29),
  [268] = {.entry = {.count = 1, .reusable = false}}, REDUCE(aux_sym__define_value_repeat1, 1),
  [270] = {.entry = {.count = 1, .reusable = true}}, REDUCE(aux_sym__define_value_repeat1, 1),
  [272] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_tuple_repeat1, 2), SHIFT_REPEAT(35),
  [275] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_tuple_repeat1, 2), SHIFT_REPEAT(3),
  [278] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_tuple_repeat1, 2), SHIFT_REPEAT(114),
  [281] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_tuple_repeat1, 2), SHIFT_REPEAT(35),
  [284] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_tuple_repeat1, 2), SHIFT_REPEAT(259),
  [287] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_tuple_repeat1, 2), SHIFT_REPEAT(299),
  [290] = {.entry = {.count = 1, .reusable = true}}, REDUCE(aux_sym_tuple_repeat1, 2),
  [292] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_range, 4),
  [294] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_range, 4),
  [296] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_block, 2),
  [298] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_block, 2),
  [300] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_tuple, 3),
  [302] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_tuple, 3),
  [304] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_collision, 10),
  [306] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_collision, 10),
  [308] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_excitation, 10),
  [310] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_excitation, 10),
  [312] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym__string_literal_single, 3),
  [314] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym__string_literal_single, 3),
  [316] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_range, 6),
  [318] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_range, 6),
  [320] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_string_literal, 1),
  [322] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_string_literal, 1),
  [324] = {.entry = {.count = 1, .reusable = true}}, SHIFT(35),
  [326] = {.entry = {.count = 1, .reusable = true}}, SHIFT(3),
  [328] = {.entry = {.count = 1, .reusable = false}}, SHIFT(35),
  [330] = {.entry = {.count = 1, .reusable = true}}, SHIFT(127),
  [332] = {.entry = {.count = 1, .reusable = true}}, SHIFT(128),
  [334] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_collision, 16),
  [336] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_collision, 16),
  [338] = {.entry = {.count = 1, .reusable = true}}, SHIFT(88),
  [340] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_range, 5),
  [342] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_range, 5),
  [344] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_list, 3),
  [346] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_list, 3),
  [348] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym__string_literal_double, 3),
  [350] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym__string_literal_double, 3),
  [352] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_new, 4),
  [354] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_new, 4),
  [356] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_block, 3),
  [358] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_block, 3),
  [360] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_unit, 1),
  [362] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_unit, 1),
  [364] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_dimension, 2),
  [366] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_dimension, 2),
  [368] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_collision, 7),
  [370] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_collision, 7),
  [372] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_assignment, 3),
  [374] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_assignment, 3),
  [376] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_boolean, 1),
  [378] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_boolean, 1),
  [380] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_generate, 3),
  [382] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_generate, 3),
  [384] = {.entry = {.count = 1, .reusable = true}}, SHIFT(49),
  [386] = {.entry = {.count = 1, .reusable = true}}, SHIFT(38),
  [388] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_rate, 1),
  [390] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_rate, 1),
  [392] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_new, 3),
  [394] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_new, 3),
  [396] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_associative_new, 6),
  [398] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_associative_new, 6),
  [400] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_arrhenius, 5),
  [402] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_arrhenius, 5),
  [404] = {.entry = {.count = 1, .reusable = true}}, SHIFT(84),
  [406] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_generate, 4),
  [408] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_generate, 4),
  [410] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_janev, 11),
  [412] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_janev, 11),
  [414] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_reaction, 6),
  [416] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_reaction, 6),
  [418] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_associative_new, 5),
  [420] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_associative_new, 5),
  [422] = {.entry = {.count = 1, .reusable = false}}, SHIFT(59),
  [424] = {.entry = {.count = 1, .reusable = true}}, SHIFT(47),
  [426] = {.entry = {.count = 1, .reusable = false}}, SHIFT(47),
  [428] = {.entry = {.count = 1, .reusable = true}}, SHIFT(65),
  [430] = {.entry = {.count = 1, .reusable = false}}, SHIFT(65),
  [432] = {.entry = {.count = 1, .reusable = true}}, SHIFT(45),
  [434] = {.entry = {.count = 1, .reusable = false}}, SHIFT(45),
  [436] = {.entry = {.count = 1, .reusable = true}}, SHIFT(60),
  [438] = {.entry = {.count = 1, .reusable = false}}, SHIFT(60),
  [440] = {.entry = {.count = 1, .reusable = true}}, SHIFT(44),
  [442] = {.entry = {.count = 1, .reusable = false}}, SHIFT(44),
  [444] = {.entry = {.count = 1, .reusable = true}}, REDUCE(aux_sym_obj_key_repeat1, 1),
  [446] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_define, 3),
  [448] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_define, 3),
  [450] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_include, 2),
  [452] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_include, 2),
  [454] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_ifxdef, 4),
  [456] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_ifxdef, 4),
  [458] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_ifxdef, 5),
  [460] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_ifxdef, 5),
  [462] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_ifxdef, 3),
  [464] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_ifxdef, 3),
  [466] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym__define_value, 2),
  [468] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym__define_value, 2),
  [470] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_get, 2),
  [472] = {.entry = {.count = 1, .reusable = false}}, REDUCE(sym_get, 2),
  [474] = {.entry = {.count = 1, .reusable = false}}, SHIFT(140),
  [476] = {.entry = {.count = 1, .reusable = true}}, SHIFT(149),
  [478] = {.entry = {.count = 1, .reusable = false}}, SHIFT(280),
  [480] = {.entry = {.count = 1, .reusable = false}}, SHIFT(279),
  [482] = {.entry = {.count = 1, .reusable = false}}, SHIFT(275),
  [484] = {.entry = {.count = 1, .reusable = false}}, SHIFT(137),
  [486] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_obj_key_repeat1, 2), SHIFT_REPEAT(139),
  [489] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_obj_key_repeat1, 2), SHIFT_REPEAT(144),
  [492] = {.entry = {.count = 2, .reusable = false}}, REDUCE(aux_sym_obj_key_repeat1, 2), SHIFT_REPEAT(140),
  [495] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_obj_key_repeat1, 2), SHIFT_REPEAT(149),
  [498] = {.entry = {.count = 1, .reusable = true}}, SHIFT(139),
  [500] = {.entry = {.count = 1, .reusable = false}}, SHIFT(312),
  [502] = {.entry = {.count = 1, .reusable = false}}, SHIFT(318),
  [504] = {.entry = {.count = 1, .reusable = false}}, SHIFT(331),
  [506] = {.entry = {.count = 1, .reusable = true}}, SHIFT(161),
  [508] = {.entry = {.count = 1, .reusable = true}}, SHIFT(292),
  [510] = {.entry = {.count = 1, .reusable = true}}, SHIFT(293),
  [512] = {.entry = {.count = 1, .reusable = true}}, SHIFT(26),
  [514] = {.entry = {.count = 1, .reusable = true}}, REDUCE(aux_sym_chems_repeat1, 2),
  [516] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_chems_repeat1, 2), SHIFT_REPEAT(176),
  [519] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_chems_repeat1, 2), SHIFT_REPEAT(224),
  [522] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_chems, 2),
  [524] = {.entry = {.count = 1, .reusable = true}}, SHIFT(176),
  [526] = {.entry = {.count = 1, .reusable = true}}, SHIFT(224),
  [528] = {.entry = {.count = 1, .reusable = true}}, SHIFT(158),
  [530] = {.entry = {.count = 1, .reusable = true}}, SHIFT(23),
  [532] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_chems, 1),
  [534] = {.entry = {.count = 1, .reusable = true}}, SHIFT(141),
  [536] = {.entry = {.count = 1, .reusable = true}}, SHIFT(22),
  [538] = {.entry = {.count = 1, .reusable = true}}, SHIFT(96),
  [540] = {.entry = {.count = 1, .reusable = true}}, SHIFT(122),
  [542] = {.entry = {.count = 1, .reusable = true}}, SHIFT(238),
  [544] = {.entry = {.count = 1, .reusable = true}}, SHIFT(250),
  [546] = {.entry = {.count = 1, .reusable = true}}, SHIFT(135),
  [548] = {.entry = {.count = 1, .reusable = true}}, SHIFT(246),
  [550] = {.entry = {.count = 1, .reusable = true}}, SHIFT(248),
  [552] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym__pterm, 2),
  [554] = {.entry = {.count = 1, .reusable = true}}, SHIFT(310),
  [556] = {.entry = {.count = 1, .reusable = true}}, SHIFT(337),
  [558] = {.entry = {.count = 1, .reusable = true}}, SHIFT(285),
  [560] = {.entry = {.count = 1, .reusable = true}}, SHIFT(296),
  [562] = {.entry = {.count = 1, .reusable = true}}, SHIFT(304),
  [564] = {.entry = {.count = 1, .reusable = true}}, SHIFT(335),
  [566] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym__nterm, 2),
  [568] = {.entry = {.count = 1, .reusable = true}}, SHIFT(267),
  [570] = {.entry = {.count = 1, .reusable = true}}, SHIFT(48),
  [572] = {.entry = {.count = 1, .reusable = true}}, SHIFT(251),
  [574] = {.entry = {.count = 1, .reusable = true}}, SHIFT(106),
  [576] = {.entry = {.count = 1, .reusable = true}}, SHIFT(287),
  [578] = {.entry = {.count = 1, .reusable = true}}, SHIFT(172),
  [580] = {.entry = {.count = 1, .reusable = true}}, SHIFT(257),
  [582] = {.entry = {.count = 1, .reusable = true}}, SHIFT(36),
  [584] = {.entry = {.count = 1, .reusable = true}}, SHIFT(148),
  [586] = {.entry = {.count = 1, .reusable = true}}, SHIFT(320),
  [588] = {.entry = {.count = 1, .reusable = true}}, SHIFT(170),
  [590] = {.entry = {.count = 1, .reusable = true}}, SHIFT(163),
  [592] = {.entry = {.count = 1, .reusable = true}}, SHIFT(329),
  [594] = {.entry = {.count = 1, .reusable = true}}, SHIFT(93),
  [596] = {.entry = {.count = 1, .reusable = true}}, SHIFT(216),
  [598] = {.entry = {.count = 1, .reusable = true}}, SHIFT(162),
  [600] = {.entry = {.count = 1, .reusable = true}}, SHIFT(186),
  [602] = {.entry = {.count = 1, .reusable = true}}, SHIFT(222),
  [604] = {.entry = {.count = 1, .reusable = true}}, SHIFT(334),
  [606] = {.entry = {.count = 1, .reusable = true}}, SHIFT(181),
  [608] = {.entry = {.count = 1, .reusable = true}}, SHIFT(192),
  [610] = {.entry = {.count = 1, .reusable = true}}, SHIFT(183),
  [612] = {.entry = {.count = 1, .reusable = true}}, SHIFT(194),
  [614] = {.entry = {.count = 1, .reusable = true}}, SHIFT(325),
  [616] = {.entry = {.count = 1, .reusable = true}}, SHIFT(185),
  [618] = {.entry = {.count = 1, .reusable = true}}, SHIFT(330),
  [620] = {.entry = {.count = 1, .reusable = true}}, SHIFT(91),
  [622] = {.entry = {.count = 1, .reusable = true}}, SHIFT(81),
  [624] = {.entry = {.count = 1, .reusable = true}}, SHIFT(196),
  [626] = {.entry = {.count = 1, .reusable = true}}, SHIFT(188),
  [628] = {.entry = {.count = 1, .reusable = true}}, SHIFT(198),
  [630] = {.entry = {.count = 1, .reusable = true}}, SHIFT(167),
  [632] = {.entry = {.count = 1, .reusable = true}}, SHIFT(199),
  [634] = {.entry = {.count = 1, .reusable = true}}, SHIFT(111),
  [636] = {.entry = {.count = 1, .reusable = true}}, SHIFT(201),
  [638] = {.entry = {.count = 1, .reusable = true}}, SHIFT(204),
  [640] = {.entry = {.count = 1, .reusable = true}}, SHIFT(193),
  [642] = {.entry = {.count = 1, .reusable = true}}, SHIFT(205),
  [644] = {.entry = {.count = 1, .reusable = true}}, SHIFT(326),
  [646] = {.entry = {.count = 1, .reusable = true}}, SHIFT(309),
  [648] = {.entry = {.count = 1, .reusable = true}}, SHIFT(206),
  [650] = {.entry = {.count = 1, .reusable = true}}, SHIFT(178),
  [652] = {.entry = {.count = 1, .reusable = true}}, SHIFT(212),
  [654] = {.entry = {.count = 1, .reusable = true}}, SHIFT(197),
  [656] = {.entry = {.count = 1, .reusable = true}}, SHIFT(200),
  [658] = {.entry = {.count = 1, .reusable = true}}, SHIFT(202),
  [660] = {.entry = {.count = 1, .reusable = true}}, REDUCE(aux_sym_full_formula_repeat1, 2),
  [662] = {.entry = {.count = 2, .reusable = true}}, REDUCE(aux_sym_full_formula_repeat1, 2), SHIFT_REPEAT(175),
  [665] = {.entry = {.count = 1, .reusable = true}}, SHIFT(217),
  [667] = {.entry = {.count = 1, .reusable = true}}, SHIFT(207),
  [669] = {.entry = {.count = 1, .reusable = true}}, SHIFT(164),
  [671] = {.entry = {.count = 1, .reusable = true}}, SHIFT(228),
  [673] = {.entry = {.count = 1, .reusable = true}}, SHIFT(113),
  [675] = {.entry = {.count = 1, .reusable = true}}, SHIFT(311),
  [677] = {.entry = {.count = 1, .reusable = true}}, SHIFT(315),
  [679] = {.entry = {.count = 1, .reusable = true}}, SHIFT(218),
  [681] = {.entry = {.count = 1, .reusable = true}}, SHIFT(225),
  [683] = {.entry = {.count = 1, .reusable = true}}, SHIFT(67),
  [685] = {.entry = {.count = 1, .reusable = true}}, SHIFT(208),
  [687] = {.entry = {.count = 1, .reusable = true}}, SHIFT(209),
  [689] = {.entry = {.count = 1, .reusable = true}}, SHIFT(239),
  [691] = {.entry = {.count = 1, .reusable = true}}, SHIFT(175),
  [693] = {.entry = {.count = 1, .reusable = true}}, SHIFT(289),
  [695] = {.entry = {.count = 1, .reusable = true}}, SHIFT(168),
  [697] = {.entry = {.count = 1, .reusable = true}}, SHIFT(64),
  [699] = {.entry = {.count = 1, .reusable = true}}, SHIFT(159),
  [701] = {.entry = {.count = 1, .reusable = true}}, SHIFT(39),
  [703] = {.entry = {.count = 1, .reusable = true}}, SHIFT(191),
  [705] = {.entry = {.count = 1, .reusable = true}}, SHIFT(55),
  [707] = {.entry = {.count = 1, .reusable = true}}, SHIFT(306),
  [709] = {.entry = {.count = 1, .reusable = true}}, SHIFT(336),
  [711] = {.entry = {.count = 1, .reusable = true}}, SHIFT(100),
  [713] = {.entry = {.count = 1, .reusable = true}}, SHIFT(46),
  [715] = {.entry = {.count = 1, .reusable = true}}, SHIFT(245),
  [717] = {.entry = {.count = 1, .reusable = true}}, SHIFT(173),
  [719] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_full_formula, 4),
  [721] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_sub_formula, 3),
  [723] = {.entry = {.count = 1, .reusable = true}}, REDUCE(sym_full_formula, 3),
  [725] = {.entry = {.count = 1, .reusable = true}}, SHIFT(305),
  [727] = {.entry = {.count = 1, .reusable = true}}, SHIFT(266),
  [729] = {.entry = {.count = 1, .reusable = true}},  ACCEPT_INPUT(),
  [731] = {.entry = {.count = 1, .reusable = true}}, SHIFT(104),
  [733] = {.entry = {.count = 1, .reusable = true}}, SHIFT(117),
  [735] = {.entry = {.count = 1, .reusable = true}}, SHIFT(215),
  [737] = {.entry = {.count = 1, .reusable = true}}, SHIFT(230),
  [739] = {.entry = {.count = 1, .reusable = true}}, SHIFT(297),
  [741] = {.entry = {.count = 1, .reusable = true}}, SHIFT(268),
  [743] = {.entry = {.count = 1, .reusable = true}}, SHIFT(298),
  [745] = {.entry = {.count = 1, .reusable = true}}, SHIFT(235),
  [747] = {.entry = {.count = 1, .reusable = true}}, SHIFT(237),
  [749] = {.entry = {.count = 1, .reusable = true}}, SHIFT(290),
  [751] = {.entry = {.count = 1, .reusable = true}}, SHIFT(78),
  [753] = {.entry = {.count = 1, .reusable = true}}, SHIFT(174),
  [755] = {.entry = {.count = 1, .reusable = true}}, SHIFT(42),
  [757] = {.entry = {.count = 1, .reusable = true}}, SHIFT(95),
  [759] = {.entry = {.count = 1, .reusable = true}}, SHIFT(94),
  [761] = {.entry = {.count = 1, .reusable = true}}, SHIFT(136),
  [763] = {.entry = {.count = 1, .reusable = true}}, SHIFT(130),
  [765] = {.entry = {.count = 1, .reusable = true}}, SHIFT(6),
  [767] = {.entry = {.count = 1, .reusable = true}}, SHIFT(13),
  [769] = {.entry = {.count = 1, .reusable = true}}, SHIFT(79),
  [771] = {.entry = {.count = 1, .reusable = true}}, SHIFT(260),
  [773] = {.entry = {.count = 1, .reusable = true}}, SHIFT(31),
  [775] = {.entry = {.count = 1, .reusable = true}}, SHIFT(286),
  [777] = {.entry = {.count = 1, .reusable = true}}, SHIFT(317),
  [779] = {.entry = {.count = 1, .reusable = true}}, SHIFT(263),
  [781] = {.entry = {.count = 1, .reusable = true}}, SHIFT(243),
  [783] = {.entry = {.count = 1, .reusable = true}}, SHIFT(265),
  [785] = {.entry = {.count = 1, .reusable = true}}, SHIFT(274),
  [787] = {.entry = {.count = 1, .reusable = true}}, SHIFT(30),
  [789] = {.entry = {.count = 1, .reusable = true}}, SHIFT(182),
  [791] = {.entry = {.count = 1, .reusable = true}}, SHIFT(220),
  [793] = {.entry = {.count = 1, .reusable = true}}, SHIFT(41),
  [795] = {.entry = {.count = 1, .reusable = true}}, SHIFT(169),
  [797] = {.entry = {.count = 1, .reusable = true}}, SHIFT(227),
  [799] = {.entry = {.count = 1, .reusable = true}}, SHIFT(255),
  [801] = {.entry = {.count = 1, .reusable = true}}, SHIFT(231),
  [803] = {.entry = {.count = 1, .reusable = true}}, SHIFT(261),
  [805] = {.entry = {.count = 1, .reusable = true}}, SHIFT(262),
  [807] = {.entry = {.count = 1, .reusable = true}}, SHIFT(16),
  [809] = {.entry = {.count = 1, .reusable = true}}, SHIFT(219),
  [811] = {.entry = {.count = 1, .reusable = true}}, SHIFT(316),
  [813] = {.entry = {.count = 1, .reusable = true}}, SHIFT(258),
  [815] = {.entry = {.count = 1, .reusable = true}}, SHIFT(269),
  [817] = {.entry = {.count = 1, .reusable = true}}, SHIFT(270),
  [819] = {.entry = {.count = 1, .reusable = true}}, SHIFT(271),
  [821] = {.entry = {.count = 1, .reusable = true}}, SHIFT(332),
  [823] = {.entry = {.count = 1, .reusable = true}}, SHIFT(214),
  [825] = {.entry = {.count = 1, .reusable = true}}, SHIFT(171),
  [827] = {.entry = {.count = 1, .reusable = true}}, SHIFT(288),
  [829] = {.entry = {.count = 1, .reusable = true}}, SHIFT(213),
  [831] = {.entry = {.count = 1, .reusable = true}}, SHIFT(223),
  [833] = {.entry = {.count = 1, .reusable = true}}, SHIFT(291),
  [835] = {.entry = {.count = 1, .reusable = true}}, SHIFT(211),
  [837] = {.entry = {.count = 1, .reusable = true}}, SHIFT(322),
  [839] = {.entry = {.count = 1, .reusable = true}}, SHIFT(236),
  [841] = {.entry = {.count = 1, .reusable = true}}, SHIFT(303),
  [843] = {.entry = {.count = 1, .reusable = true}}, SHIFT(307),
  [845] = {.entry = {.count = 1, .reusable = true}}, SHIFT(324),
  [847] = {.entry = {.count = 1, .reusable = true}}, SHIFT(195),
  [849] = {.entry = {.count = 1, .reusable = true}}, SHIFT(50),
  [851] = {.entry = {.count = 1, .reusable = true}}, SHIFT(313),
  [853] = {.entry = {.count = 1, .reusable = true}}, SHIFT(314),
  [855] = {.entry = {.count = 1, .reusable = true}}, SHIFT(177),
  [857] = {.entry = {.count = 1, .reusable = true}}, SHIFT(203),
  [859] = {.entry = {.count = 1, .reusable = true}}, SHIFT(190),
  [861] = {.entry = {.count = 1, .reusable = true}}, SHIFT(319),
  [863] = {.entry = {.count = 1, .reusable = true}}, SHIFT(189),
  [865] = {.entry = {.count = 1, .reusable = true}}, SHIFT(321),
  [867] = {.entry = {.count = 1, .reusable = true}}, SHIFT(323),
  [869] = {.entry = {.count = 1, .reusable = true}}, SHIFT(7),
  [871] = {.entry = {.count = 1, .reusable = true}}, SHIFT(86),
  [873] = {.entry = {.count = 1, .reusable = true}}, SHIFT(284),
  [875] = {.entry = {.count = 1, .reusable = true}}, SHIFT(327),
  [877] = {.entry = {.count = 1, .reusable = true}}, SHIFT(229),
  [879] = {.entry = {.count = 1, .reusable = true}}, SHIFT(187),
  [881] = {.entry = {.count = 1, .reusable = true}}, SHIFT(256),
  [883] = {.entry = {.count = 1, .reusable = true}}, SHIFT(184),
  [885] = {.entry = {.count = 1, .reusable = true}}, SHIFT(333),
  [887] = {.entry = {.count = 1, .reusable = true}}, SHIFT(180),
};

#ifdef __cplusplus
extern "C" {
#endif
#ifdef _WIN32
#define extern __declspec(dllexport)
#endif

extern const TSLanguage *tree_sitter_turbowave(void) {
  static const TSLanguage language = {
    .version = LANGUAGE_VERSION,
    .symbol_count = SYMBOL_COUNT,
    .alias_count = ALIAS_COUNT,
    .token_count = TOKEN_COUNT,
    .external_token_count = EXTERNAL_TOKEN_COUNT,
    .state_count = STATE_COUNT,
    .large_state_count = LARGE_STATE_COUNT,
    .production_id_count = PRODUCTION_ID_COUNT,
    .field_count = FIELD_COUNT,
    .max_alias_sequence_length = MAX_ALIAS_SEQUENCE_LENGTH,
    .parse_table = &ts_parse_table[0][0],
    .small_parse_table = ts_small_parse_table,
    .small_parse_table_map = ts_small_parse_table_map,
    .parse_actions = ts_parse_actions,
    .symbol_names = ts_symbol_names,
    .symbol_metadata = ts_symbol_metadata,
    .public_symbol_map = ts_symbol_map,
    .alias_map = ts_non_terminal_alias_map,
    .alias_sequences = &ts_alias_sequences[0][0],
    .lex_modes = ts_lex_modes,
    .lex_fn = ts_lex,
  };
  return &language;
}
#ifdef __cplusplus
}
#endif