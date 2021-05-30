#include "ParentVesselReconstruction.h"
// qt
#include <QFile>

// vtk
#include "vtkCleanPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataNormals.h"
#include "vtkDelaunay3D.h"
#include "vtkSplineFilter.h"

// vmtk
#include "vtkvmtkCapPolyData.h"
#include "vtkvmtkInternalTetrahedraExtractor.h"
#include "vtkvmtkVoronoiDiagram3D.h"
#include "vtkvmtkSimplifyVoronoiDiagram.h"

#include "vtkXMLPolyDataReader.h"


// me
#include "interactorStylePickCenterline.h"

ParentVesselReconstruction::ParentVesselReconstruction(QObject* parent)
{
	
}

ParentVesselReconstruction::~ParentVesselReconstruction()
{

}

void ParentVesselReconstruction::SetSource(vtkPolyData *source)
{
	m_source->DeepCopy(source);
}

void ParentVesselReconstruction::SetCenterline(vtkPolyData *centerline)
{
	m_centerline->DeepCopy(centerline);
}

vtkPolyData * ParentVesselReconstruction::GetSource()
{
	return m_source;
}

vtkPolyData * ParentVesselReconstruction::GetClippedCenterline()
{
	return m_clippedCenterline;
}

vtkPolyData * ParentVesselReconstruction::GetClippedVoronoiDiagram()
{
	return m_clippedVoronoiDiagram;
}

vtkPolyData * ParentVesselReconstruction::GetVoronoiDiagram()
{
	return m_voronoiDiagram;
}

vtkPolyData * ParentVesselReconstruction::GetCenterline()
{
	return m_centerline;
}

void ParentVesselReconstruction::Run()
{
	std::cout << "capping input surface..." << std::endl;
	// capping
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(m_source);
	cleaner->Update();

	vtkSmartPointer<vtkTriangleFilter> triangulator = vtkSmartPointer<vtkTriangleFilter>::New();
	triangulator->SetInputData(cleaner->GetOutput());
	triangulator->PassLinesOff();
	triangulator->PassVertsOff();
	triangulator->Update();

	m_source->DeepCopy(triangulator->GetOutput());

	// quick load voronoi diagram for debug
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName("Z:\\data\\intracranial\\data_ESASIS_followup\\medical\\055\\baseline\\clipped_voronoi.vtp");
	reader->SetFileName("D:/Projects/Parent-Vessel-Reconstruction/Data/voronoi.vtp");
	reader->Update();
	m_voronoiDiagram->DeepCopy(reader->GetOutput());

	//this->ComputeVoronoiDiagram();

	std::cout << "subdivide centerline..." << std::endl;
	vtkSmartPointer<vtkSplineFilter> splineFilter = vtkSmartPointer<vtkSplineFilter>::New();
	splineFilter->SetInputData(m_centerline);
	splineFilter->SetSubdivideToLength();
	splineFilter->SetLength(0.1);
	splineFilter->Update();

	m_centerline->DeepCopy(splineFilter->GetOutput());

	this->SeedPicker();
}

void ParentVesselReconstruction::SeedPicker()
{
	vtkSmartPointer<vtkPolyDataMapper> surfaceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	surfaceMapper->SetInputData(m_source);
	surfaceMapper->SetScalarVisibility(false);
	vtkSmartPointer<vtkActor> surfaceActor = vtkSmartPointer<vtkActor>::New();
	surfaceActor->SetMapper(surfaceMapper);
	surfaceActor->GetProperty()->SetOpacity(0.3);

	vtkSmartPointer<vtkPolyDataMapper> centerlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	centerlineMapper->SetInputData(m_centerline);
	centerlineMapper->SetScalarVisibility(false);
	vtkSmartPointer<vtkActor> centerlineActor = vtkSmartPointer<vtkActor>::New();
	centerlineActor->SetMapper(centerlineMapper);

	// Create a sphere
	m_sphere->SetCenter(m_centerline->GetPoint(0));
	m_sphere->SetRadius(0.5);

	vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	sphereMapper->SetInputConnection(m_sphere->GetOutputPort());

	vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
	sphereActor->SetMapper(sphereMapper);
	sphereActor->GetProperty()->SetColor(1, 0, 0);

	// put the actor into render window
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	ren->AddActor(surfaceActor);
	ren->AddActor(centerlineActor);
	ren->AddActor(sphereActor);

	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(ren);
	renWin->SetSize(1024, 768);

	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<MouseInteractorStylePickCenterline> style = vtkSmartPointer<MouseInteractorStylePickCenterline>::New();
	iren->SetRenderWindow(renWin);
	iren->SetInteractorStyle(style);
	style->SetCenterline(m_centerline);
	style->SetClippedCenterline(m_clippedCenterline);
	style->SetSphere(m_sphere);
	style->SetVoronoiDiagram(m_voronoiDiagram);
	style->SetClippedVoronoiDiagram(m_clippedVoronoiDiagram);
	style->SetReconstructedSurface(m_source);

	iren->Initialize();
	renWin->Render();
	ren->ResetCamera();
	iren->Start();
}

void ParentVesselReconstruction::ComputeVoronoiDiagram()
{
	std::cout << "start to compute Voronoi diagram..." << std::endl;

	vtkSmartPointer<vtkvmtkCapPolyData> capper = vtkSmartPointer<vtkvmtkCapPolyData>::New();
	capper->SetInputData(m_source);
	capper->Update();

	vtkSmartPointer<vtkPolyDataNormals> surfaceNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
	surfaceNormals->SetInputData(capper->GetOutput());
	surfaceNormals->SplittingOff();
	surfaceNormals->AutoOrientNormalsOn();
	surfaceNormals->SetFlipNormals(false);
	surfaceNormals->ComputePointNormalsOn();
	surfaceNormals->ConsistencyOn();
	surfaceNormals->Update();

	std::cout << "performing delaunay tesselation..." << std::endl;

	vtkSmartPointer<vtkUnstructuredGrid> delaunayTessellation = vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkDelaunay3D> delaunayTessellator = vtkSmartPointer<vtkDelaunay3D>::New();
	delaunayTessellator->CreateDefaultLocator();
	delaunayTessellator->SetInputConnection(surfaceNormals->GetOutputPort());
	delaunayTessellator->SetTolerance(0.001);
	delaunayTessellator->Update();
	delaunayTessellation->DeepCopy(delaunayTessellator->GetOutput());

	vtkDataArray* normalsArray = surfaceNormals->GetOutput()->GetPointData()->GetNormals();
	delaunayTessellation->GetPointData()->AddArray(normalsArray);

	std::cout << "extracting internal tetrahedra..." << std::endl;

	vtkSmartPointer<vtkvmtkInternalTetrahedraExtractor> internalTetrahedraExtractor = vtkSmartPointer<vtkvmtkInternalTetrahedraExtractor>::New();
	internalTetrahedraExtractor->SetInputData(delaunayTessellation);
	internalTetrahedraExtractor->SetOutwardNormalsArrayName(normalsArray->GetName());
	internalTetrahedraExtractor->RemoveSubresolutionTetrahedraOn();
	internalTetrahedraExtractor->SetSubresolutionFactor(1.0);
	internalTetrahedraExtractor->SetSurface(surfaceNormals->GetOutput());

	if (capper->GetCapCenterIds()->GetNumberOfIds() > 0)
	{
		internalTetrahedraExtractor->UseCapsOn();
		internalTetrahedraExtractor->SetCapCenterIds(capper->GetCapCenterIds());
		internalTetrahedraExtractor->Update();
	}

	delaunayTessellation->DeepCopy(internalTetrahedraExtractor->GetOutput());

	std::cout << "computing Voronoi diagram..." << std::endl;

	vtkSmartPointer<vtkvmtkVoronoiDiagram3D> voronoiDiagramFilter = vtkSmartPointer<vtkvmtkVoronoiDiagram3D>::New();
	voronoiDiagramFilter->SetInputData(delaunayTessellation);
	voronoiDiagramFilter->SetRadiusArrayName("Radius");
	voronoiDiagramFilter->Update();

	std::cout << "simplifying Voronoi diagram..." << std::endl;

	vtkSmartPointer<vtkvmtkSimplifyVoronoiDiagram> voronoiDiagramSimplifier = vtkSmartPointer<vtkvmtkSimplifyVoronoiDiagram>::New();
	voronoiDiagramSimplifier->SetInputData(voronoiDiagramFilter->GetOutput());
	voronoiDiagramSimplifier->SetUnremovablePointIds(voronoiDiagramFilter->GetPoleIds());
	voronoiDiagramSimplifier->Update();

	m_voronoiDiagram->DeepCopy(voronoiDiagramSimplifier->GetOutput());
}