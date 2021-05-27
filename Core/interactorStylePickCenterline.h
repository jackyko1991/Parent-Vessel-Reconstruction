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
	void SetCenterline(vtkPolyData* centerline);
	void SetSphere(vtkSphereSource* sphere);
	void SetVoronoiDiagram(vtkPolyData* voronoiDiagram);
	vtkPolyData* GetClippedCenterline();
	vtkPolyData* GetOutputCenterline();
	vtkPolyData* GetVoronoiDiagram();

private:
	vtkSmartPointer<vtkSphereSource> m_sphere = NULL;
	vtkPolyData* m_centerline = NULL;
	vtkSmartPointer<vtkPolyData> m_clippedCenterline = NULL;
	vtkSmartPointer<vtkPolyData> m_outputCenterline = NULL;
	vtkSmartPointer<vtkPolyData> m_voronoiDiagram = NULL;
	vtkSmartPointer<vtkPolyData> m_clippedVoronoiDiagram = NULL;
	vtkSmartPointer<vtkPolyData> m_outputVoronoiDiagram = NULL;
	vtkSmartPointer<vtkPolyData> m_normalized_centerline = NULL;
	vtkSmartPointer<vtkKdTreePointLocator> m_kDTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	
	int m_pickedPointId = 0;
	double m_clipDistance = 2.5;
	double m_comThreshold = 1.5;
	QList< QPair<vtkActor* ,vtkActor* >> m_clipPlaneActorList;
	bool m_smoothVoronoiDiagram = false;

	void ClipPlaneUpdate();
	void CreateClipPlaneActor(vtkUnstructuredGrid*, int, vtkActor*);
	void ClipCenterline();
	void ClipVoronoiDiagram();
	void InterpoldateVoronoiDiagram();
	void ExtractCylindricInterpolationVoronoiDiagram(int, int, double, vtkPolyData*, vtkPolyData*, vtkPolyData*);

};

#endif