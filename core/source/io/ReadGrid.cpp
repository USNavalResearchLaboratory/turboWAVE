module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"

export module read_grid;
import base;
import metric_space;
import driver;
import logger;
import input;

/// This class reads the grid block and some outer directives from the input file.
/// This data is then used to help initialize Task and MetricSpace, after which the
/// GridReader can be thrown out.
/// Neither MetricSpace nor Task read anything directly themselves.
export class GridReader
{
	tw::input::DirectiveReader outerDirectives;
	tw::bc::par bc0[4],bc1[4];
	tw::vec4 min_spacing,max_spacing,critical_spacing;
	tw::vec4 maxWindowPosition,solutionVelocity;
	bool neutralize,movingWindow;
	
	tw::input::DirectiveReader directives;
	tw::grid::geometry geo;
	tw::Int steps,req_dim[4],req_dom[4];
	tw::vec4 relativeRef,absoluteRef,spacing;
	bool adaptiveTimestep,adaptiveGrid,found;

	public:
	GridReader(tw::units native,tw::Float unitDensityCGS) {
		min_spacing = tw::vec4(tw::small_pos);
		max_spacing = tw::vec4(tw::big_pos);
		critical_spacing = tw::vec4(tw::small_pos);
		solutionVelocity = tw::vec4(1.0,0.0,0.0,0.0);
		maxWindowPosition = tw::big_pos;
		neutralize = false;
		movingWindow = false;
		bc0[1] = tw::bc::par::periodic;
		bc1[1] = tw::bc::par::periodic;
		bc0[2] = tw::bc::par::periodic;
		bc1[2] = tw::bc::par::periodic;
		bc0[3] = tw::bc::par::absorbing;
		bc1[3] = tw::bc::par::absorbing;

		outerDirectives.Add("xboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[1],&bc1[1]));
		outerDirectives.Add("yboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[2],&bc1[2]));
		outerDirectives.Add("zboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[3],&bc1[3]));
		outerDirectives.Add("dtmin",new tw::input::Float(&min_spacing[0]),false);
		outerDirectives.Add("dtmax",new tw::input::Float(&max_spacing[0]),false);
		outerDirectives.Add("dtcrit",new tw::input::Float(&critical_spacing[0]),false);
		outerDirectives.Add("maxtime",new tw::input::Float(&maxWindowPosition[0]),false);
		outerDirectives.Add("neutralize",new tw::input::Bool(&neutralize),false);
		outerDirectives.Add("window speed",new tw::input::Float(&solutionVelocity[3]),false);
		outerDirectives.Add("moving window",new tw::input::Bool(&movingWindow),false);
		outerDirectives.AttachUnits(native,unitDensityCGS);

		adaptiveGrid = adaptiveTimestep = found = false;
		geo = tw::grid::cartesian;
		absoluteRef = 0.0;

		directives.Add("steps",new tw::input::Int(&steps),true);
		directives.Add("origin",new tw::input::Numbers<tw::Float>(&relativeRef[0],4));
		directives.Add("shift",new tw::input::Numbers<tw::Float>(&absoluteRef[0],4),false);
		directives.Add("cell size",new tw::input::Numbers<tw::Float>(&spacing[0],4));
		directives.Add("adaptive timestep",new tw::input::Bool(&adaptiveTimestep),false);
		directives.Add("adaptive grid",new tw::input::Bool(&adaptiveGrid),false);
		directives.Add("dimensions",new tw::input::Numbers<tw::Int>(&req_dim[0],4));
		directives.Add("decomposition",new tw::input::Numbers<tw::Int>(&req_dom[0],4));
		std::map<std::string,tw::grid::geometry> geo_map = {{"cartesian",tw::grid::cartesian},{"cylindrical",tw::grid::cylindrical},{"spherical",tw::grid::spherical}};
		directives.Add("geometry",new tw::input::Enums<tw::grid::geometry>(geo_map,&geo),false);
		directives.AttachUnits(native,unitDensityCGS);
	}
	bool ReadOuter(const TSTreeCursor *curs,const std::string& src) {
		return outerDirectives.ReadNext(curs,src);
	}
	bool ReadGridBlock(const TSTreeCursor *curs0,const std::string& src) {
		auto curs = tw::input::Cursor(curs0);
		if (tw::input::node_kind(curs.get())=="new") {
			ts_tree_cursor_goto_first_child(curs.get()); // `new` token
			ts_tree_cursor_goto_next_sibling(curs.get());
			if (tw::input::node_text(curs.get(),src)=="grid") {
				ts_tree_cursor_goto_next_sibling(curs.get());
				if (tw::input::node_kind(curs.get())=="block") {
					ts_tree_cursor_goto_first_child(curs.get());
					tw::input::next_named_node(curs.get(),true);
					directives.ReadAll(curs.get(),src);
					directives.ThrowErrorIfMissingKeys("grid");
					if (req_dom[0]!=1)
						tw::input::ThrowParsingError(curs.get(),src,"Time decomposition must be 1");
					found = true;
					return true;
				} else {
					tw::input::ThrowParsingError(curs.get(),src,"Ill-formed grid block");
				}
			}
		}
		return false;
	}
	void UpdateTask(Task& tsk) {
		// these are subject to transformation before being fed back into Task::Initialize by Simulation
		tsk.periodic[1] = bc0[1]==tw::bc::par::periodic ? 1 : 0;
		tsk.periodic[2] = bc0[2]==tw::bc::par::periodic ? 1 : 0;
		tsk.periodic[3] = bc0[3]==tw::bc::par::periodic ? 1 : 0;
		for (auto i=0;i<4;i++) {
			tsk.domains[i] = req_dom[i];
		}
	}
	void UpdateSpace(MetricSpace& ms) {
		// whatever not handled by Resize method.
		// n.b. this is hard-coded for time evolutions.
		for (auto i=0; i<4; i++) {
			ms.bc0[i] = bc0[i];
			ms.bc1[i] = bc1[i];
		}
		ms.adaptiveGrid = adaptiveGrid;
		ms.adaptiveTimestep = adaptiveTimestep;
		ms.neutralize = neutralize;
		if (movingWindow && !outerDirectives.TestKey("window speed")) {
			solutionVelocity[3] = 1.0;
		}
		ms.SetupInternalEvolution(0,spacing[0],steps,0);
		ms.ChangeStepSizeControls(0,min_spacing[0],max_spacing[0],critical_spacing[0]);
		ms.SetupBoundaryEvolution(solutionVelocity,maxWindowPosition);
	}
	tw::Int Steps() {
		return steps;
	}
	tw::vec4 GlobalCorner() {
		return absoluteRef - relativeRef*GlobalSize();
	}
	tw::vec4 GlobalSize() {
		return tw::vec4(
			spacing[0] * req_dim[0],
			spacing[1] * req_dim[1],
			spacing[2] * req_dim[2],
			spacing[3] * req_dim[3]
		);
	}
	tw::idx4 GlobalDims() {
		return tw::idx4(req_dim[0],req_dim[1],req_dim[2],req_dim[3]);
	}
	tw::grid::geometry Geometry() {
		return geo;
	}
	bool FoundGrid() {
		return found;
	}
};
