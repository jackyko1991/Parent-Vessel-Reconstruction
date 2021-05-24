#include "ParentVesselReconstruction.h"
// qt
#include <QFile>

// vtk
#include "vtkSTLReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPolyDataReader.h"
#include "observe_error.h""
#include "vtkCleanPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

// vmtk
#include "vtkvmtkCapPolyData.h"

// me
#include "interactorStyleCenterline.h"

ParentVesselReconstruction::ParentVesselReconstruction(QObject* parent)
{
}

ParentVesselReconstruction::~ParentVesselReconstruction()
{

}

void ParentVesselReconstruction::SetSourceFilePath(QString path)
{
	m_sourceFile.setFile(path);
}

void ParentVesselReconstruction::SetOutputFilePath(QString path)
{
	m_outputFile.setFile(path);
}

void ParentVesselReconstruction::Run()
{
	// check input file existence
	if (!m_sourceFile.exists())
	{
		std::cout << "Source file not exist!!!" << std::endl;
		return;
	}

	// read surface
	vtkSmartPointer<ErrorObserver> errorObserver = vtkSmartPointer<ErrorObserver>::New();

	if (m_sourceFile.suffix() == "vtp" || m_sourceFile.suffix() == "VTP")
	{
		vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
		reader->SetFileName(m_sourceFile.absoluteFilePath().toStdString().c_str());
		reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
		reader->Update();
		m_source->DeepCopy(reader->GetOutput());
	}
	else if (m_sourceFile.suffix() == "stl" || m_sourceFile.suffix() == "STL")
	{
		vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(m_sourceFile.absoluteFilePath().toStdString().c_str());
		reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
		reader->Update();
		m_source->DeepCopy(reader->GetOutput());

	}
	else if (m_sourceFile.suffix() == "vtk" || m_sourceFile.suffix() == "VTK")
	{
		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
		reader->SetFileName(m_sourceFile.absoluteFilePath().toStdString().c_str());
		reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
		m_source->DeepCopy(reader->GetOutput());
	}
	else
	{
		std::cerr << "Invalid input data type, only accept STL, VTP or VTK files" << std::endl;
		return;
	}

	// capping
	vtkSmartPointer<vtkvmtkCapPolyData> capper = vtkSmartPointer<vtkvmtkCapPolyData>::New();
	capper->SetInputData(m_source);
	capper->Update();

	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(capper->GetOutput());
	cleaner->Update();

	vtkSmartPointer<vtkTriangleFilter> triangulator = vtkSmartPointer<vtkTriangleFilter>::New();
	triangulator->SetInputData(cleaner->GetOutput());
	triangulator->PassLinesOff();
	triangulator->PassVertsOff();
	triangulator->Update();

	m_source->DeepCopy(triangulator->GetOutput());

	this->SeedPicker();
}

void ParentVesselReconstruction::SeedPicker()
{
	vtkSmartPointer<vtkPolyDataMapper> surfaceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	surfaceMapper->SetInputData(m_source);
	vtkSmartPointer<vtkActor> surfaceActor = vtkSmartPointer<vtkActor>::New();
	surfaceActor->SetMapper(surfaceMapper);

	vtkSmartPointer<vtkPolyDataMapper> centerlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	centerlineMapper->SetInputData(m_centerline);
	vtkSmartPointer<vtkActor> centerlineActor = vtkSmartPointer<vtkActor>::New();
	centerlineActor->SetMapper(centerlineMapper);

	// put the actor into render window
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	ren->AddActor(surfaceActor);
	ren->AddActor(centerlineActor);

	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(ren);

	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<MouseInteractorStyleCenterline> style = vtkSmartPointer<MouseInteractorStyleCenterline>::New();
	iren->SetInteractorStyle(style);
	iren->SetRenderWindow(renWin);
	style->SetSurface(m_source);
	style->SetCenterline(m_centerline);

	iren->Initialize();
	renWin->Render();
	ren->ResetCamera();
	iren->Start();
}
