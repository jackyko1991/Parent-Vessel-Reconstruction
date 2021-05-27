#ifndef __PARENT_VESSEL_RECONSTRUCTION_H__
#define __PARENT_VESSEL_RECONSTRUCTION_H__

// qt
#include <QObject>
#include <QFileInfo>

// vtk
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <vtkSphereSource.h>

class ParentVesselReconstruction : public QObject
{
	Q_OBJECT

public:
	explicit ParentVesselReconstruction(QObject* parent = 0);
	~ParentVesselReconstruction();

	void SetSourceFilePath(QString);
	void SetCenterlineFilePath(QString);
	void SetOutputFilePath(QString);
	void Run();

private:
	void SeedPicker();
	void ComputeVoronoiDiagram();
	
	QFileInfo m_sourceFile;
	QFileInfo m_centerlineFile;
	QFileInfo m_outputFile;

	vtkSmartPointer<vtkPolyData> m_source = vtkSmartPointer<vtkPolyData>::New() ;
	vtkSmartPointer<vtkPolyData> m_centerline = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> m_clipped_centerline = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> m_voronoiDiagram = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkSphereSource > m_sphere = vtkSmartPointer<vtkSphereSource>::New();
};

#endif