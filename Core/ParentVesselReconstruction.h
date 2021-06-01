#ifndef __PARENT_VESSEL_RECONSTRUCTION_H__
#define __PARENT_VESSEL_RECONSTRUCTION_H__

// qt
#include <QObject>
#include <QFileInfo>

// vtk
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <vtkSphereSource.h>
#include "vtkPolyData.h"

class ParentVesselReconstruction : public QObject
{
	Q_OBJECT

public:
	explicit ParentVesselReconstruction(QObject* parent = 0);
	~ParentVesselReconstruction();

	void SetSource(vtkPolyData*);
	void SetCenterline(vtkPolyData*);
	void SetVoronoiSmooth(double);
	void SetCenterOfMassThreshold(double);
	vtkPolyData* GetSource();
	vtkPolyData* GetClippedCenterline();
	vtkPolyData* GetClippedVoronoiDiagram();
	vtkPolyData* GetVoronoiDiagram();
	vtkPolyData* GetCenterline();

	void Run();

private:
	void SeedPicker();
	void ComputeVoronoiDiagram();

	double m_voronoiSmooth = 0.6;
	double m_comThreshold = 1.5;
	vtkSmartPointer<vtkPolyData> m_source = vtkSmartPointer<vtkPolyData>::New() ;
	vtkSmartPointer<vtkPolyData> m_centerline = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> m_clippedCenterline = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> m_voronoiDiagram = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> m_clippedVoronoiDiagram = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkSphereSource > m_sphere = vtkSmartPointer<vtkSphereSource>::New();
};

#endif