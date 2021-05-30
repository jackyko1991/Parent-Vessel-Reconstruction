#ifndef __INTERACTOR_STYLE_PICK_CENTERLINE_H__
#define __INTERACTOR_STYLE_PICK_CENTERLINE_H__

#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPointData.h>
#include <vtkObjectFactory.h>
#include <vtkPointPicker.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkCoordinate.h>
#include <vtkKdTreePointLocator.h>
#include <vtkSphereSource.h>
#include <vtkRenderer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkActorCollection.h>

#include <Centerline.h>

#include <QMap>
#include <QList>
#include <QPair>

class MouseInteractorStylePickCenterline : public vtkInteractorStyleTrackballCamera
{
public:
	static MouseInteractorStylePickCenterline* New();
	vtkTypeMacro(MouseInteractorStylePickCenterline, vtkInteractorStyleTrackballCamera);

	void OnKeyPress();
	void SetReconstructedSurface(vtkPolyData*);
	void SetCenterline(vtkPolyData* centerline);
	void SetClippedCenterline(vtkPolyData* clippedCenterline);
	void SetSphere(vtkSphereSource* sphere);
	void SetVoronoiDiagram(vtkPolyData* voronoiDiagram);
	void SetClippedVoronoiDiagram(vtkPolyData* clippedVoronoiDiagram);

private:
	vtkSmartPointer<vtkSphereSource> m_sphere = NULL;
	vtkSmartPointer<vtkPolyData> m_centerline = NULL;
	vtkSmartPointer<vtkPolyData> m_clippedCenterline = NULL;
	vtkSmartPointer<vtkPolyData> m_voronoiDiagram = NULL;
	vtkSmartPointer<vtkPolyData> m_clippedVoronoiDiagram = NULL;
	vtkSmartPointer<vtkPolyData> m_normalized_centerline = NULL;
	vtkSmartPointer<vtkPolyData> m_reconstructedSurface = NULL;
	vtkSmartPointer<vtkKdTreePointLocator> m_kDTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	
	int m_pickedPointId = 0;
	double m_clipDistance = 2.5;
	double m_comThreshold = 1.5;
	double m_smoothFactor = 0.4;
	QList< QPair<vtkActor* ,vtkActor* >> m_clipPlaneActorList;
	bool m_smoothVoronoiDiagram = false;
	double m_pointCloudDensity = 200;
	int m_clippedCenterlinesCount = 1;

	void ClipPlaneUpdate();
	void CreateClipPlaneActor(vtkUnstructuredGrid*, int, vtkActor*);
	void ClipCenterline();
	void ClipVoronoiDiagram();
	void InterpoldateVoronoiDiagram();
	vtkPolyData* ExtractCylindricInterpolationVoronoiDiagram(vtkPolyData*);
};

#endif