#include "meta_base.h"

MetricSpace::MetricSpace()
{
	geo = tw::grid::cartesian;
	car = 1.0;
	cyl = sph = 0.0;
	I3x3[0][0] = 1;
	I3x3[0][1] = 0;
	I3x3[0][2] = 0;
	I3x3[1][0] = 0;
	I3x3[1][1] = 1;
	I3x3[1][2] = 0;
	I3x3[2][0] = 0;
	I3x3[2][1] = 0;
	I3x3[2][2] = 1;
	units = tw::UnitConverter(tw::units::plasma,1e19);
	#ifdef USE_OPENCL
	metricsBuffer = NULL;
	#endif
}

MetricSpace::~MetricSpace()
{
	#ifdef USE_OPENCL
	if (metricsBuffer!=NULL)
	{
		clReleaseMemObject(metricsBuffer);
		clReleaseMemObject(stripBuffer[0]);
		clReleaseMemObject(stripBuffer[1]);
		clReleaseMemObject(stripBuffer[2]);
		clReleaseMemObject(stripBuffer[3]);
	}
	#endif
}

void MetricSpace::AttachUnits(tw::units sys,tw::Float unitDensityCGS)
{
	units = tw::UnitConverter(sys,unitDensityCGS);
}

#ifdef USE_OPENCL

void MetricSpace::InitializeMetricsBuffer(cl_context ctx,tw::Float dt)
{
	cl_int err;
	tw::Float metrics[12] = { dt , spacing.x , spacing.y , spacing.z , 0.0 , corner.x , corner.y , corner.z , car , cyl , sph , 0.0 };
	tw::Int tStrip[7] = { 0,0,0,0,0,0 }; // not generally used; can signal a special operation
	tw::Int xStrip[7] = { 1,1,0,0,encodingStride[1],dim[1] }; // which-axis , x-axis? , y-axis? , z-axis? , stride , dim (there is deliberate redundancy)
	tw::Int yStrip[7] = { 2,0,1,0,encodingStride[2],dim[2] };
	tw::Int zStrip[7] = { 3,0,0,1,encodingStride[3],dim[3] };

	metricsBuffer = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Float)*12,metrics,&err);
	stripBuffer[0] = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Int)*6,tStrip,&err);
	stripBuffer[1] = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Int)*6,xStrip,&err);
	stripBuffer[2] = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Int)*6,yStrip,&err);
	stripBuffer[3] = clCreateBuffer(ctx,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(tw::Int)*6,zStrip,&err);

	if (err)
		throw tw::FatalError("Device buffer failure in MetricSpace.");
}

void MetricSpace::StripUpdateProtocol(cl_kernel k,cl_command_queue q,tw::Int axis,tw::Int stripArgument)
{
	size_t global_work_range[3] = { (size_t)num[1],(size_t)num[2],(size_t)num[3] };
	global_work_range[axis-1] = 1;
	clSetKernelArg(k,stripArgument,sizeof(stripBuffer[axis]),&stripBuffer[axis]);
	clEnqueueNDRangeKernel(q,k,3,NULL,global_work_range,NULL,0,NULL,NULL);
	clFinish(q);
}

#endif

void MetricSpace::Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size,tw::Int ghostCellLayers)
{
	DiscreteSpace::Resize(x,y,z,corner,size,ghostCellLayers);

	// Metric arrays have ghost cell layers even when dim=1.

	mlb[1] = 1 - layers[0];
	mlb[2] = 1 - layers[0];
	mlb[3] = 1 - layers[0];

	mub[1] = dim[1] + layers[0];
	mub[2] = dim[2] + layers[0];
	mub[3] = dim[3] + layers[0];

	mnum[1] = dim[1]+2*layers[0];
	mnum[2] = dim[2]+2*layers[0];
	mnum[3] = dim[3]+2*layers[0];

	mnum[0] = mnum[1];
	if (mnum[2]>mnum[0])
		mnum[0] = mnum[2];
	if (mnum[3]>mnum[0])
		mnum[0] = mnum[3];

	gpos.resize(4*mnum[0]); // node coordinates
	width.resize(4*mnum[0]); // parametric cell size (not an arc length in general)

	cell_area_x.resize(mnum[1]*8);
	cell_arc_x.resize(mnum[1]*6);
	cell_area_z.resize(mnum[3]*8);
	cell_arc_z.resize(mnum[3]*6);

	// Setup up a uniform grid in parameter space
	// Non-uniform grids can be externally imposed after this by overwriting the cell widths

	for (tw::Int i=mlb[1];i<=mub[1];i++)
		dX(i,1) = spacing.x;
		for (tw::Int i=mlb[2];i<=mub[2];i++)
		dX(i,2) = spacing.y;
		for (tw::Int i=mlb[3];i<=mub[3];i++)
		dX(i,3) = spacing.z;
}

void MetricSpace::Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size)
{
	Resize(x,y,z,corner,size,2);  // default to 2 ghost cell layers
}

void MetricSpace::SetupPositionArrays()
{
	// Assumes corner and width have already been setup.
	// Typically not called directly, instead call Set*Geometry()
	tw::Int i,ax;

	for (ax=1;ax<=3;ax++)
	{
		X(mlb[ax],ax) = corner[ax-1] - dX(mlb[ax]+1,ax) - 0.5*dX(mlb[ax],ax);
		for (i=mlb[ax]+1;i<=mub[ax];i++)
			X(i,ax) = X(i-1,ax) + 0.5*dX(i-1,ax) + 0.5*dX(i,ax);
	}
}

void MetricSpace::SetCartesianGeometry()
{
	SetupPositionArrays();
	geo = tw::grid::cartesian;
	car = 1.0; cyl = 0.0; sph = 0.0;
	tw::Int i,sx=mnum[1],sz=mnum[3];
	tw::Float *xpos = &gpos[mnum[0]*1];
	tw::Float *ypos = &gpos[mnum[0]*2];
	tw::Float *zpos = &gpos[mnum[0]*3];
	tw::Float *xwidth = &width[mnum[0]*1];
	tw::Float *ywidth = &width[mnum[0]*2];
	tw::Float *zwidth = &width[mnum[0]*3];
	for (i=0;i<sx;i++)
	{
		cell_area_x[i + 0*sx] = xwidth[i]*ywidth[2];
		cell_area_x[i + 1*sx] = ywidth[2];
		cell_area_x[i + 2*sx] = xwidth[i];
		cell_area_x[i + 3*sx] = xwidth[i]*ywidth[2];
		cell_arc_x[i + 1*sx] = ywidth[2];
		cell_arc_x[i + 2*sx] = 1.0;
	}
	for (i=0;i<sz;i++)
	{
		cell_area_z[i + 0*sz] = zwidth[i];
		cell_area_z[i + 1*sz] = zwidth[i];
		cell_area_z[i + 2*sz] = zwidth[i];
		cell_area_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 0*sz] = 1.0;
		cell_arc_z[i + 1*sz] = 1.0;
	}
	for (i=1;i<sx;i++)
	{
		cell_arc_x[i + 0*sx] = xpos[i] - xpos[i-1];
	}
	for (i=1;i<sz;i++)
	{
		cell_arc_z[i + 2*sz] = zpos[i] - zpos[i-1];
	}
	for (i=1;i<sx;i++)
	{
		cell_area_x[i + 4*sx] = (xpos[i]-xpos[i-1])*ywidth[2];
		cell_area_x[i + 6*sx] = (xpos[i]-xpos[i-1]);
		cell_area_x[i + 7*sx] = (xpos[i]-xpos[i-1])*ywidth[2];
	}
	for (i=1;i<sz;i++)
	{
		cell_area_z[i + 4*sz] = zpos[i]-zpos[i-1];
		cell_area_z[i + 5*sz] = zpos[i]-zpos[i-1];
		cell_area_z[i + 6*sz] = zpos[i]-zpos[i-1];
	}
	for (i=0;i<sx;i++)
	{
		cell_area_x[i + 5*sx] = ywidth[2];
		cell_arc_x[i + 3*sx] = xwidth[i];
		cell_arc_x[i + 4*sx] = ywidth[2];
		cell_arc_x[i + 5*sx] = 1.0;
	}
	for (i=0;i<sz;i++)
	{
		cell_area_z[i + 7*sz] = 1.0;
		cell_arc_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 4*sz] = 1.0;
		cell_arc_z[i + 5*sz] = zwidth[i];
	}
}

void MetricSpace::SetCylindricalGeometry()
{
	// we are using x = r, y = phi, z = z
	SetupPositionArrays();
	geo = tw::grid::cylindrical;
	car = 0.0; cyl = 1.0; sph = 0.0;
	tw::Int i,sx=mnum[1],sz=mnum[3];
	tw::Float *xpos = &gpos[mnum[0]*1];
	tw::Float *ypos = &gpos[mnum[0]*2];
	tw::Float *zpos = &gpos[mnum[0]*3];
	tw::Float *xwidth = &width[mnum[0]*1];
	tw::Float *ywidth = &width[mnum[0]*2];
	tw::Float *zwidth = &width[mnum[0]*3];
	for (i=0;i<sx;i++)
	{
		const tw::Float x0 = xpos[i]; const tw::Float x1 = xpos[i] - 0.5*xwidth[i];
		cell_area_x[i + 0*sx] = x0*xwidth[i]*ywidth[2];
		cell_area_x[i + 1*sx] = x1*ywidth[2];
		cell_area_x[i + 2*sx] = xwidth[i];
		cell_area_x[i + 3*sx] = x0*xwidth[i]*ywidth[2];
		cell_arc_x[i + 1*sx] = x0*ywidth[2];
		cell_arc_x[i + 2*sx] = 1.0;
	}
	for (i=0;i<sz;i++)
	{
		cell_area_z[i + 0*sz] = zwidth[i];
		cell_area_z[i + 1*sz] = zwidth[i];
		cell_area_z[i + 2*sz] = zwidth[i];
		cell_area_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 0*sz] = 1.0;
		cell_arc_z[i + 1*sz] = 1.0;
	}
	for (i=1;i<sx;i++)
	{
		cell_arc_x[i + 0*sx] = xpos[i] - xpos[i-1];
	}
	for (i=1;i<sz;i++)
	{
		cell_arc_z[i + 2*sz] = zpos[i] - zpos[i-1];
	}
	for (i=1;i<sx;i++)
	{
		const tw::Float x0 = xpos[i] - 0.5*xwidth[i];
		cell_area_x[i + 4*sx] = x0*(xpos[i]-xpos[i-1])*ywidth[2];
		cell_area_x[i + 6*sx] = (xpos[i]-xpos[i-1]);
		cell_area_x[i + 7*sx] = x0*(xpos[i]-xpos[i-1])*ywidth[2];
	}
	for (i=1;i<sz;i++)
	{
		cell_area_z[i + 4*sz] = zpos[i]-zpos[i-1];
		cell_area_z[i + 5*sz] = zpos[i]-zpos[i-1];
		cell_area_z[i + 6*sz] = zpos[i]-zpos[i-1];
	}
	for (i=0;i<sx;i++)
	{
		const tw::Float x0 = xpos[i] - 0.5*xwidth[i]; const tw::Float x1 = xpos[i];
		cell_area_x[i + 5*sx] = x1*ywidth[2];
		cell_arc_x[i + 3*sx] = xwidth[i];
		cell_arc_x[i + 4*sx] = x0*ywidth[2];
		cell_arc_x[i + 5*sx] = 1.0;
	}
	for (i=0;i<sz;i++)
	{
		cell_area_z[i + 7*sz] = 1.0;
		cell_arc_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 4*sz] = 1.0;
		cell_arc_z[i + 5*sz] = zwidth[i];
	}

	// Don't allow any wall area to strictly collapse
	for (i=0;i<sx;i++)
		for (tw::Int c=0;c<8;c++)
			if (cell_area_x[i + c*sx]==0.0)
				cell_area_x[i + c*sx] = tw::small_pos;
}

void MetricSpace::SetSphericalGeometry()
{
	// we are using q1 = rho, q2 = phi (atan(y/x)), q3 = theta (acos(z/rho))
	SetupPositionArrays();
	geo = tw::grid::spherical;
	car = 0.0; cyl = 0.0; sph = 1.0;
	tw::Int i,sx=mnum[1],sz=mnum[3];
	tw::Float *xpos = &gpos[mnum[0]*1];
	tw::Float *ypos = &gpos[mnum[0]*2];
	tw::Float *zpos = &gpos[mnum[0]*3];
	tw::Float *xwidth = &width[mnum[0]*1];
	tw::Float *ywidth = &width[mnum[0]*2];
	tw::Float *zwidth = &width[mnum[0]*3];
	for (i=0;i<sx;i++)
	{
		const tw::Float x0 = xpos[i];
		const tw::Float dx = xwidth[i]; const tw::Float dy = ywidth[2];
		const tw::Float x1 = x0 - 0.5*dx;
		cell_area_x[i + 0*sx] = dx*dy*(dx*dx/6.0 + 2.0*x0*x0);
		cell_area_x[i + 1*sx] = 2.0*x1*x1*dy;
		cell_area_x[i + 2*sx] = dx*x0;
		cell_area_x[i + 3*sx] = dx*dy*x0;
		cell_arc_x[i + 1*sx] = x0*ywidth[2];
		cell_arc_x[i + 2*sx] = x0;
	}
	for (i=0;i<sz;i++)
	{
		const tw::Float z0 = zpos[i];
		const tw::Float dz = zwidth[i];
		const tw::Float z1 = z0 - 0.5*dz;
		cell_area_z[i + 0*sz] = sin(0.5*dz)*sin(z0);
		cell_area_z[i + 1*sz] = sin(0.5*dz)*sin(z0);
		cell_area_z[i + 2*sz] = dz;
		cell_area_z[i + 3*sz] = sin(z1);
		cell_arc_z[i + 0*sz] = 1.0;
		cell_arc_z[i + 1*sz] = sin(z0);
	}
	for (i=1;i<sx;i++)
	{
		cell_arc_x[i + 0*sx] = xpos[i] - xpos[i-1];
	}
	for (i=1;i<sz;i++)
	{
		cell_arc_z[i + 2*sz] = zpos[i] - zpos[i-1];
	}
	for (i=1;i<sx;i++)
	{
		const tw::Float x0 = xpos[i] - 0.5*xwidth[i];
		const tw::Float dx = xpos[i] - xpos[i-1]; const tw::Float dy = ywidth[2];
		cell_area_x[i + 4*sx] = dx*dy*(dx*dx/6.0 + 2.0*x0*x0);
		cell_area_x[i + 6*sx] = dx*x0;
		cell_area_x[i + 7*sx] = dx*dy*x0;
	}
	for (i=1;i<sz;i++)
	{
		const tw::Float z0 = zpos[i] - 0.5*zwidth[i];
		const tw::Float dz = zpos[i] - zpos[i-1];
		cell_area_z[i + 4*sz] = sin(0.5*dz)*sin(z0);
		cell_area_z[i + 5*sz] = sin(0.5*dz)*sin(z0);
		cell_area_z[i + 6*sz] = dz;
	}
	for (i=0;i<sx;i++)
	{
		const tw::Float x0 = xpos[i] - 0.5*xwidth[i];
 		const tw::Float dy = ywidth[2];
		const tw::Float x1 = xpos[i];
		cell_area_x[i + 5*sx] = 2.0*x1*x1*dy;
		cell_arc_x[i + 3*sx] = xwidth[i];
		cell_arc_x[i + 4*sx] = x0*dy;
		cell_arc_x[i + 5*sx] = x0;
	}
	for (i=0;i<sz;i++)
	{
		const tw::Float z0 = zpos[i] - 0.5*zwidth[i];
		const tw::Float z1 = zpos[i];
		cell_area_z[i + 7*sz] = sin(z1);
		cell_arc_z[i + 3*sz] = 1.0;
		cell_arc_z[i + 4*sz] = sin(z0);
		cell_arc_z[i + 5*sz] = zwidth[i];
	}

	// Don't allow any wall area to strictly collapse
	// N.b. however tw::small_pos^2 = underflow
	for (i=0;i<sx;i++)
		for (tw::Int c=0;c<8;c++)
			if (cell_area_x[i + c*sx]==0.0)
				cell_area_x[i + c*sx] = tw::small_pos;
	for (i=0;i<sz;i++)
		for (tw::Int c=0;c<8;c++)
			if (cell_area_z[i + c*sz]==0.0)
				cell_area_z[i + c*sz] = tw::small_pos;
}

void MetricSpace::ReadCheckpoint(std::ifstream& inFile)
{
	DiscreteSpace::ReadCheckpoint(inFile);
	Resize(dim[1],dim[2],dim[3],corner,size,layers[0]);
	inFile.read((char *)mnum,sizeof(mnum));
	inFile.read((char *)mlb,sizeof(mlb));
	inFile.read((char *)mub,sizeof(mub));
	inFile.read((char *)&car,sizeof(tw::Float));
	inFile.read((char *)&cyl,sizeof(tw::Float));
	inFile.read((char *)&sph,sizeof(tw::Float));
	inFile.read((char *)&gpos[0],sizeof(tw::Float)*gpos.size());
	inFile.read((char *)&width[0],sizeof(tw::Float)*width.size());
	inFile.read((char *)&cell_area_x[0],sizeof(tw::Float)*cell_area_x.size());
	inFile.read((char *)&cell_arc_x[0],sizeof(tw::Float)*cell_arc_x.size());
	inFile.read((char *)&cell_area_z[0],sizeof(tw::Float)*cell_area_z.size());
	inFile.read((char *)&cell_arc_z[0],sizeof(tw::Float)*cell_arc_z.size());
}

void MetricSpace::WriteCheckpoint(std::ofstream& outFile)
{
	DiscreteSpace::WriteCheckpoint(outFile);
	outFile.write((char *)mnum,sizeof(mnum));
	outFile.write((char *)mlb,sizeof(mlb));
	outFile.write((char *)mub,sizeof(mub));
	outFile.write((char *)&car,sizeof(tw::Float));
	outFile.write((char *)&cyl,sizeof(tw::Float));
	outFile.write((char *)&sph,sizeof(tw::Float));
	outFile.write((char *)&gpos[0],sizeof(tw::Float)*gpos.size());
	outFile.write((char *)&width[0],sizeof(tw::Float)*width.size());
	outFile.write((char *)&cell_area_x[0],sizeof(tw::Float)*cell_area_x.size());
	outFile.write((char *)&cell_arc_x[0],sizeof(tw::Float)*cell_arc_x.size());
	outFile.write((char *)&cell_area_z[0],sizeof(tw::Float)*cell_area_z.size());
	outFile.write((char *)&cell_arc_z[0],sizeof(tw::Float)*cell_arc_z.size());
}
