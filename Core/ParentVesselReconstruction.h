#ifndef __PARENT_VESSEL_RECONSTRUCTION_H__
#define __PARENT_VESSEL_RECONSTRUCTION_H__

// qt
#include <QObject>
#include <QFileInfo>

// vtk
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

class ParentVesselReconstruction : public QObject
{
	Q_OBJECT

public:
	explicit ParentVesselReconstruction(QObject* parent = 0);
	~ParentVesselReconstruction();

	void SetSourceFilePath(QString);
	void SetOutputFilePath(QString);
	void Run();

private:
	void SeedPicker();
	
	QFileInfo m_sourceFile;
	QFileInfo m_outputFile;

	vtkSmartPointer<vtkPolyData> m_source = vtkSmartPointer<vtkPolyData>::New() ;
	vtkSmartPointer<vtkPolyData> m_centerline = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> m_clipped_centerline = vtkSmartPointer<vtkPolyData>::New();


};

#endif