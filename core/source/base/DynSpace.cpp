module;

#include "tw_includes.h"
#include "tw_logger.h"

export module dyn_space;
import base;
import pic_primitives;
export import tensor;
export import tasks;
import static_space;
import logger;

/// This object adds dynamical quantities to StaticSpace, but stops short of adding metric information.
/// An important task is to track the way fixed reference points move relative to StaticSpace.
export struct DynSpace : StaticSpace
{
	protected:

	tw::Int stepNow,stepsToTake;
	/// There is a marker in the inertial frame, this is its velocity.
	/// Ordinary initial value problems will have (1,0,0,0), while "moving windows" will have (1,0,0,1).
	tw::vec4 solutionVelocity;
	/// There is a marker in the inertial frame, this is its position.
	/// The marker will not in general be in perfect synchronism with the window.
	tw::vec4 solutionPosition;
	tw::vec4 altSolutionPosition;
	tw::vec4 windowPosition; ///< the actual position of the computational window
	tw::vec4 altWindowPosition;
	tw::vec4 maxWindowPosition; ///< used to determine when simulation should stop
	tw::vec4 corner; ///< position where all coordinates are minimized on the local domain
	tw::vec4 size; ///< length of the local domain along each axis
	tw::vec4 globalCorner; ///< position where all coordinates are minimized on the global domain
	tw::vec4 globalSize; ///< length of the global domain along each axis
	tw::vec4 min_spacing; ///< minimum spacing to use for adaptive grids (including time levels)
	tw::vec4 max_spacing; ///< maximum spacing to use for adaptive grids (including time levels)
	tw::vec4 critical_spacing; ///< threshold that triggers some action when we have an adaptive grid
	/// number of interior cells along the given axis, summed over all domains on the low-side of this domain
	tw::node4 lowSideCells;
	/// number of interior cells along the given axis, summed over all domains
	tw::node4 globalCells;

	public:

	/// Create an empty `DynSpace`, the normal pattern is to do this and then `Resize`.
	DynSpace();
	/// Create a `DynSpace` with a single time node and purely *local* coordinates
	DynSpace(const tw::node5& dim,const tw::vec4& corner,const tw::vec4& size,const tw::node5& packing,const tw::node4& ghostCellLayers);
	/// Change the topology and coordinates:  `gdim`, `gcorner`, and `gsize` describe global window that exists in distributed
	/// memory at any one time, as distinct from the union of all windows that occur during the system evolution.
	/// This makes MPI calls on `task` to get the domain information.
	void Resize(Task *task,const tw::node5& gdim,const tw::vec4& gcorner,const tw::vec4& gsize,const tw::node5& packing,const tw::node4& ghostCellLayers);
	void SetupTimeInfo(tw::Float dt0) { spacing[0] = dt0; freq[0] = 1.0/dt0; }
	void UpdateWindow(const DynSpace& src) {
		solutionPosition = src.solutionPosition;
		solutionVelocity = src.solutionVelocity;
		windowPosition = src.windowPosition;
		altSolutionPosition = src.altSolutionPosition;
		altWindowPosition = src.altWindowPosition;
		maxWindowPosition = src.maxWindowPosition;
	}
	void SetPrimitiveWithPosition(Primitive& q,const tw::vec4& pos) const;
	tw::vec4 PositionFromPrimitive(const Primitive& q) const;
	tw::Int GlobalDim(const tw::Int& ax) const { return globalCells[ax]; }
	tw::Int GlobalCellIndex(const tw::Int& idx,const tw::Int& ax) const { return idx + lowSideCells[ax]; }
	tw::Int LocalCellIndex(const tw::Int& idx,const tw::Int& ax) const { return idx - lowSideCells[ax]; }
	tw::Int StepNow() const { return stepNow; }
	tw::Int StepsToTake() const { return stepsToTake; }
	tw::vec4 PhysicalSize() { return size; }
	tw::vec4 GlobalPhysicalSize() { return globalSize; }
	tw::vec4 Corner() { return corner; }
	tw::vec4 GlobalCorner() { return globalCorner; }
	tw::Float SolutionPos(const tw::Int& ax) const { return solutionPosition[ax]; }
	tw::Float WindowPos(const tw::Int& ax) const { return windowPosition[ax]; }
	tw::Float MaxWindowPos(const tw::Int& ax) const { return maxWindowPosition[ax]; }
	tw::Float CriticalSpacing(const tw::Int& ax) const { return critical_spacing[ax]; }
	tw::Float MinSpacing(const tw::Int& ax) const { return min_spacing[ax]; }
	tw::Float MaxSpacing(const tw::Int& ax) const { return max_spacing[ax]; }
	bool IsPointWithinInterior(const tw::vec4& P);
	void GetWeights(weights_3D* weights,const tw::vec4& P);

	void ReadCheckpoint(std::ifstream& inFile) {
		inFile.read((char *)this,sizeof(*this));
	}
	void WriteCheckpoint(std::ofstream& outFile) {
		outFile.write((char *)this,sizeof(*this));
	}
};

DynSpace::DynSpace()
{
	min_spacing = tw::vec4(tw::small_pos);
	max_spacing = tw::vec4(tw::big_pos);
	critical_spacing = tw::vec4(tw::small_pos);
	solutionVelocity = tw::vec4(1.0,0.0,0.0,0.0);
}

DynSpace::DynSpace(const tw::node5& dim,const tw::vec4& corner,const tw::vec4& size,const tw::node5& packing,const tw::node4& ghostCellLayers):
	StaticSpace(dim,size,packing,ghostCellLayers)
{
	globalCorner = corner;
	globalSize = size;
	for (tw::Int i=0;i<4;i++) {
		globalCells[i] = dim[i];
		lowSideCells[i] = 0;
		this->size[i] = size[i];
		this->corner[i] = corner[i];
	}
}

void DynSpace::Resize(Task *task,const tw::node5& gdim,const tw::vec4& gcorner,const tw::vec4& gsize,const tw::node5& packing,const tw::node4& ghostCellLayers)
{
	tw::node4 domainIndex = task->strip[0].Get_coords4();
	tw::node5 domainCount { 1, task->strip[1].Get_size(), task->strip[2].Get_size(), task->strip[3].Get_size(), 1 };
	if (gdim[0] != 1) {
		logger::WARN(std::format("time dimension was not 1 ({})",gdim[0]));
	}
	StaticSpace::Resize(domainCount,gdim,gsize,packing,ghostCellLayers);
	logger::DEBUG(std::format("resize domain to {}x{}x{}x{}",dim[0],dim[1],dim[2],dim[3]));

	globalCorner = gcorner;
	globalSize = gsize;
	for (tw::Int i=0;i<4;i++) {
		globalCells[i] = gdim[i];
		lowSideCells[i] = domainIndex[i]*dim[i];
		size[i] = dim[i]*spacing[i];
		corner[i] = gcorner[i] + domainIndex[i]*size[i];
	}
}

inline bool DynSpace::IsPointWithinInterior(const tw::vec4& P)
{
	const tw::vec4 PLoc = P - corner;
	return PLoc[1]>=0.0 && PLoc[1]<size[1] && PLoc[2]>=0.0 && PLoc[2]<size[2] && PLoc[3]>=0.0 && PLoc[3]<size[3];
}

inline tw::vec4 DynSpace::PositionFromPrimitive(const Primitive& q) const
{
	tw::Int ijk[4];
	tw::vec4 ans;
	DecodeCell(q,ijk);
	for (tw::Int i=0;i<4;i++)
		ans[i] = corner[i] + spacing[i] * (tw::Float(q.x[i]) + tw::Float(ijk[i]) - tw::Float(0.5));
	return ans;
}

inline void DynSpace::SetPrimitiveWithPosition(Primitive& q,const tw::vec4& P) const
{
	tw::Int ijk[4];
	const tw::vec4 PLoc(P - corner);
	for (tw::Int i=0;i<4;i++)
		ijk[i] = tw::Int(std::floor(PLoc[i]*freq[i])) + 1;
	q.cell = EncodeCell(ijk[0],ijk[1],ijk[2],ijk[3]);
	for (tw::Int i=0;i<4;i++)
		q.x[i] = PLoc[i]*freq[i] - tw::Float(ijk[i]) + tw::Float(0.5);
}

inline void DynSpace::GetWeights(weights_3D* weights,const tw::vec4& P)
{
	Primitive q;
	SetPrimitiveWithPosition(q,P);
	StaticSpace::GetWeights(weights,q);
}