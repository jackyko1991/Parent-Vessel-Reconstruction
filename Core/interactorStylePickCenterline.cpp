#include "interactorStylePickCenterline.h"

#include "vtkArrayCalculator.h"
#include "vtkThreshold.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkPlaneSource.h"
#include "vtkDataArray.h"
#include "vtkClipPolyData.h"
#include "vtkAppendPolyData.h"
#include "vtkParametricSpline.h"
#include "vtkParametricFunctionSource.h"
#include "vtkCleanPolyData.h"
#include "vtkGenericCell.h"
#include "vtkCenterOfMass.h"
#include "vtkGeometryFilter.h"
#include "vtkConnectivityFilter.h"
#include "vtkSphere.h"
#include "vtkNew.h"
#include "vtkCylinder.h"
#include "vtkPlane.h"
#include "vtkDoubleArray.h"
#include "vtkImplicitBoolean.h"
#include "vtkPointLocator.h"
#include "vtkPointInterpolator.h"
#include "vtkPointSource.h"
#include "vtkMarchingCubes.h"

#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"

#include "vtkvmtkPolyBallLine.h"
#include "vtkvmtkCenterlineAttributesFilter.h"
#include "vtkvmtkPolyBallModeller.h"

void MouseInteractorStylePickCenterline::OnKeyPress()
{
	vtkRenderWindowInteractor *rwi = this->Interactor;
	std::string key = rwi->GetKeySym();

	if (key == "space")
	{
		std::cout << "press space" << std::endl;
	
		this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],
			this->Interactor->GetEventPosition()[1],
			0,  // always zero.
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];

		this->Interactor->GetPicker()->GetPickPosition(picked);

		// get nearest point on centerline
		m_pickedPointId = m_kDTree->FindClosestPoint(picked);
		m_sphere->SetCenter(m_centerline->GetPoint(m_pickedPointId));

		// normalize abscissas to picked point
		vtkSmartPointer<vtkArrayCalculator> converter = vtkSmartPointer<vtkArrayCalculator>::New();
		converter->SetInputData(m_centerline);
		converter->AddScalarArrayName("Abscissas");
		double pickedPointAbscissas = m_centerline->GetPointData()->GetArray("Abscissas")->GetTuple(m_pickedPointId)[0];
		char func[100];

		std::cout << "pickedPoint abscissas: " << pickedPointAbscissas << std::endl;

		sprintf(func, "Abscissas - %f", pickedPointAbscissas);
		converter->SetFunction(func);
		converter->SetResultArrayName("NormalizedAbscissas");
		converter->Update();

		std::cout << "conversion complete" << std::endl;

		m_normalized_centerline = (vtkPolyData*)converter->GetOutput();

		this->ClipPlaneUpdate();
	}
	else if (key == "plus" || key=="equal")
	{
		std::cout << "increase distance" << std::endl;

		if (m_clipDistance < 999)
		{
			m_clipDistance += 1;
			this->ClipPlaneUpdate();
		}
	}
	else if (key == "minus")
	{
		std::cout << "reduce distance" << std::endl;
		if (m_clipDistance - 1 > 0)
		{
			m_clipDistance -= 1;
			this->ClipPlaneUpdate();
		}
	}
	else if (key == "Return")
	{
		// remove all clip plane actors
		while (m_clipPlaneActorList.size() > 0)
		{
			this->GetCurrentRenderer()->RemoveActor(m_clipPlaneActorList.back().first);
			this->GetCurrentRenderer()->RemoveActor(m_clipPlaneActorList.back().second);
			m_clipPlaneActorList.pop_back();
		}

		this->ClipCenterline();
		std::cout << "clip centerline ok" << std::endl;

		this->ClipVoronoiDiagram();
		this->InterpoldateVoronoiDiagram();

		// Reconstructing Surface from Voronoi Diagram
		vtkNew<vtkvmtkPolyBallModeller> modeller;
		modeller->SetInputData(m_outputVoronoiDiagram);
		modeller->SetRadiusArrayName("MaximumInscribedSphereRadius");
		modeller->UsePolyBallLineOff();
		int polyBallImageSize[3] = { 90,90,90 };
		modeller->SetSampleDimensions(polyBallImageSize);
		modeller->Update();

		vtkNew<vtkMarchingCubes> marchingCube;
		marchingCube->SetInputData(modeller->GetOutput());
		marchingCube->SetValue(0, 0.0);
		marchingCube->Update();

		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetInputData(marchingCube->GetOutput());
		writer->SetFileName("D:/Projects/Parent-Vessel-Reconstruction/Data/output.vtp");
		writer->Write();
		std::cout << "save complete" << std::endl;
	}

	this->GetInteractor()->Render();
}

void MouseInteractorStylePickCenterline::SetCenterline(vtkPolyData* centerline)
{
	m_centerline = centerline;
	m_kDTree->SetDataSet(m_centerline);
	m_kDTree->BuildLocator();
}

void MouseInteractorStylePickCenterline::SetSphere(vtkSphereSource * sphere)
{
	m_sphere = sphere;
}

void MouseInteractorStylePickCenterline::SetVoronoiDiagram(vtkPolyData * voronoiDiagram)
{
	m_voronoiDiagram = voronoiDiagram;
}

void MouseInteractorStylePickCenterline::ClipPlaneUpdate()
{
	int iD = m_pickedPointId;

	// threshold on normalized abscissas
	double threshold_bound[2] = {
		m_normalized_centerline->GetPointData()->GetArray("NormalizedAbscissas")->GetTuple(0)[0],
		m_normalized_centerline->GetPointData()->GetArray("NormalizedAbscissas")->GetTuple(m_normalized_centerline->GetNumberOfPoints() - 1)[0] };

	if (-1 * m_clipDistance > threshold_bound[0])
	{
		threshold_bound[0] = -1 * m_clipDistance;
	}

	if (m_clipDistance < threshold_bound[1])
	{
		threshold_bound[1] = m_clipDistance;
	}

	m_normalized_centerline->GetPointData()->SetActiveScalars("NormalizedAbscissas");
	vtkSmartPointer<vtkClipPolyData> clipper1 = vtkSmartPointer<vtkClipPolyData>::New();
	clipper1->SetValue(threshold_bound[1]);
	clipper1->SetInsideOut(true);
	clipper1->GenerateClippedOutputOn();
	clipper1->SetInputData(m_normalized_centerline);
	clipper1->Update();

	vtkSmartPointer<vtkClipPolyData> clipper2 = vtkSmartPointer<vtkClipPolyData>::New();
	clipper2->SetValue(threshold_bound[0]);
	clipper2->SetInsideOut(false);
	clipper2->GenerateClippedOutputOn();
	clipper2->SetInputData(clipper1->GetOutput());
	clipper2->Update();

	//std::cout << clipper2->GetOutput()->GetCellData()->GetArray("CenterlineIds")->GetNumberOfTuples() << std::endl;

	// remove all previous rendered clippers
	while (m_clipPlaneActorList.size() > 0)
	{
		this->GetCurrentRenderer()->RemoveActor(m_clipPlaneActorList.back().first);
		this->GetCurrentRenderer()->RemoveActor(m_clipPlaneActorList.back().second);
		m_clipPlaneActorList.pop_back();
	}

	// loop over all centerlineids
	for (int i = 0; i < clipper1->GetOutput()->GetCellData()->GetArray("CenterlineIds")->GetNumberOfTuples(); i++)
	{
		// threshold to get independent lines
		vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
		threshold->ThresholdBetween(i, i);
		threshold->SetInputData(clipper2->GetOutput());
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "CenterlineIds");
		threshold->Update();

		if (threshold->GetOutput()->GetNumberOfPoints() == 0)
			continue;

		// check centroid of the thresheld section is close to picked point
		vtkSmartPointer<vtkCenterOfMass> comFilter = vtkSmartPointer<vtkCenterOfMass>::New();
		comFilter->SetInputData(threshold->GetOutput());
		comFilter->Update();
		
		double centroid_dist =
			vtkMath::Distance2BetweenPoints(
				m_centerline->GetPoint(m_pickedPointId),
				comFilter->GetCenter()
				);
		
		if (centroid_dist > m_comThreshold)
			continue;

		// create proximal and distal clip plane pairs
		// proximal plane
		vtkSmartPointer<vtkActor> proximal_actor = vtkSmartPointer<vtkActor>::New();
		this->CreateClipPlaneActor(threshold->GetOutput(), 0, proximal_actor);
		this->GetCurrentRenderer()->AddActor(proximal_actor);

		// distal plane
		vtkSmartPointer<vtkActor> distal_actor = vtkSmartPointer<vtkActor>::New();
		this->CreateClipPlaneActor(threshold->GetOutput(), threshold->GetOutput()->GetNumberOfPoints() - 1, distal_actor);
		this->GetCurrentRenderer()->AddActor(distal_actor);

		m_clipPlaneActorList.append(QPair<vtkActor*, vtkActor*>(proximal_actor, distal_actor));
	}
}

void MouseInteractorStylePickCenterline::CreateClipPlaneActor(vtkUnstructuredGrid *centerline, int id, vtkActor *actor)
{
	double center2origin[3];
	center2origin[0] = -1 * (centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id)[0] + centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id)[0]);
	center2origin[1] = -1 * (centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id)[1] + centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id)[1]);
	center2origin[2] = -1 * (centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id)[2] + centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id)[2]);

	double center2origin_norm = vtkMath::Norm(center2origin);
	double normal_norm = vtkMath::Norm(centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id));
	double binormal_norm = vtkMath::Norm(centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id));

	double radius = centerline->GetPointData()->GetArray("Radius")->GetTuple(id)[0];
	double origin[3];

	double size_factor = 3.;

	origin[0] = centerline->GetPoint(id)[0] + size_factor * radius*center2origin[0] / center2origin_norm;
	origin[1] = centerline->GetPoint(id)[1] + size_factor * radius*center2origin[1] / center2origin_norm;
	origin[2] = centerline->GetPoint(id)[2] + size_factor * radius*center2origin[2] / center2origin_norm;

	double point1[3];
	point1[0] = origin[0] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id)[0] / normal_norm;
	point1[1] = origin[1] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id)[1] / normal_norm;
	point1[2] = origin[2] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id)[2] / normal_norm;

	double point2[3];
	point2[0] = origin[0] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id)[0] / binormal_norm;
	point2[1] = origin[1] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id)[1] / binormal_norm;
	point2[2] = origin[2] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id)[2] / binormal_norm;

	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
	plane->SetOrigin(origin);
	plane->SetPoint1(point1);
	plane->SetPoint2(point2);
	plane->Update();

	vtkSmartPointer<vtkPolyDataMapper> priximal_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	priximal_mapper->SetInputData(plane->GetOutput());

	actor->SetMapper(priximal_mapper);
	actor->GetProperty()->SetColor(1, 0, 0);
}

void MouseInteractorStylePickCenterline::ClipCenterline()
{
	// create new clipped centerline
	if (m_clippedCenterline != NULL)
	{
		m_clippedCenterline->Delete();
	}
	if (m_clippedCenterline != NULL)
	{
		m_outputCenterline->Delete();
	}

	m_clippedCenterline = vtkSmartPointer<vtkPolyData>::New();
	m_outputCenterline = vtkSmartPointer<vtkPolyData>::New();

	// threshold on normalized abscissas
	double threshold_bound[2] = {
		m_normalized_centerline->GetPointData()->GetArray("NormalizedAbscissas")->GetTuple(0)[0],
		m_normalized_centerline->GetPointData()->GetArray("NormalizedAbscissas")->GetTuple(m_normalized_centerline->GetNumberOfPoints() - 1)[0] };

	if (-1 * m_clipDistance > threshold_bound[0])
	{
		threshold_bound[0] = -1 * m_clipDistance;
	}

	if (m_clipDistance < threshold_bound[1])
	{
		threshold_bound[1] = m_clipDistance;
	}

	m_normalized_centerline->GetPointData()->SetActiveScalars("NormalizedAbscissas");

	// loop over all centerlineids
	for (int i = 0; i < m_normalized_centerline->GetCellData()->GetArray("CenterlineIds")->GetNumberOfTuples(); i++)
	{
		// threshold to get independent lines
		vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
		threshold->ThresholdBetween(i, i);
		threshold->SetInputData(m_normalized_centerline);
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "CenterlineIds");
		threshold->Update();

		// convert threshold output to vtkpolydata
		vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter> ::New();
		geometryFilter->SetInputData(threshold->GetOutput());
		geometryFilter->Update();

		if (threshold->GetOutput()->GetNumberOfPoints() == 0)
			continue;

		vtkSmartPointer<vtkPolyData> singleCenterline = geometryFilter->GetOutput();
		singleCenterline->GetPointData()->SetActiveScalars("NormalizedAbscissas");

		// compute section to clip and determine if the centroid is close enought
		vtkSmartPointer<vtkClipPolyData> clipper1 = vtkSmartPointer<vtkClipPolyData>::New();
		clipper1->SetValue(threshold_bound[1]);
		clipper1->SetInsideOut(true);
		clipper1->GenerateClippedOutputOn();
		clipper1->SetInputData(singleCenterline);
		clipper1->Update();

		vtkSmartPointer<vtkClipPolyData> clipper2 = vtkSmartPointer<vtkClipPolyData>::New();
		clipper2->SetValue(threshold_bound[0]);
		clipper2->SetInsideOut(false);
		clipper2->GenerateClippedOutputOn();
		clipper2->SetInputData(clipper1->GetOutput());
		clipper2->Update();

		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner->SetInputData(clipper2->GetOutput());
		cleaner->Update();

		vtkSmartPointer<vtkCenterOfMass> comFilter = vtkSmartPointer<vtkCenterOfMass>::New();
		comFilter->SetInputData(cleaner->GetOutput());
		comFilter->Update();

		double centroid_dist =
			vtkMath::Distance2BetweenPoints(
				m_centerline->GetPoint(m_pickedPointId),
				comFilter->GetCenter()
			);

		vtkSmartPointer<vtkAppendPolyData> appendFilter2 = vtkSmartPointer<vtkAppendPolyData>::New();
		appendFilter2->AddInputData(m_outputCenterline);

		vtkSmartPointer<vtkAppendPolyData> appendFilterClipped = vtkSmartPointer<vtkAppendPolyData>::New();
		appendFilterClipped->AddInputData(m_clippedCenterline);

		if (centroid_dist > m_comThreshold)
		{
			// direct append to output
			appendFilterClipped->AddInputData(singleCenterline);
			appendFilter2->AddInputData(singleCenterline);
		}
		else
		{
			// clip the centerline and append to output
			vtkSmartPointer<vtkClipPolyData> clipper3 = vtkSmartPointer<vtkClipPolyData>::New();
			clipper3->SetValue(threshold_bound[0]);
			clipper3->SetInsideOut(true);
			clipper3->GenerateClippedOutputOn();
			clipper3->SetInputData(singleCenterline);
			clipper3->Update();

			vtkSmartPointer<vtkClipPolyData> clipper4 = vtkSmartPointer<vtkClipPolyData>::New();
			clipper4->SetValue(threshold_bound[1]);
			clipper4->SetInsideOut(false);
			clipper4->GenerateClippedOutputOn();
			clipper4->SetInputData(singleCenterline);
			clipper4->Update();

			vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
			appendFilter->AddInputData(clipper3->GetOutput());
			appendFilter->AddInputData(clipper4->GetOutput());
			appendFilter->Update();

			vtkSmartPointer<vtkCleanPolyData> cleaner2 = vtkSmartPointer<vtkCleanPolyData>::New();
			cleaner2->SetInputData(appendFilter->GetOutput());
			cleaner2->Update();

			appendFilterClipped->AddInputData(cleaner2->GetOutput());

			// interpolation with spline
			vtkNew<vtkParametricSpline> spline;
			spline->SetPoints(cleaner2->GetOutput()->GetPoints());

			vtkNew<vtkParametricFunctionSource> functionSource;
			functionSource->SetParametricFunction(spline);
			functionSource->Update();

			appendFilter2->AddInputData(functionSource->GetOutput());
		}
		appendFilterClipped->Update();
		appendFilter2->Update();
		m_clippedCenterline->DeepCopy(appendFilterClipped->GetOutput());
		m_outputCenterline->DeepCopy(appendFilter2->GetOutput());
	}

	// recompute centerline attributes
	vtkSmartPointer<vtkvmtkCenterlineAttributesFilter> attributeFilter = vtkSmartPointer<vtkvmtkCenterlineAttributesFilter>::New();
	attributeFilter->SetInputData(m_outputCenterline);
	attributeFilter->SetAbscissasArrayName("Abscissas");
	attributeFilter->SetParallelTransportNormalsArrayName("ParallelTransportNormals");
	attributeFilter->Update();

	m_outputCenterline->DeepCopy(attributeFilter->GetOutput());
	m_centerline->DeepCopy(m_outputCenterline);
}

void MouseInteractorStylePickCenterline::ClipVoronoiDiagram()
{
	if (m_voronoiDiagram == NULL)
		return;

	std::cout << "clipping voronoi diagram" << std::endl;

	// compute connectivity on clipped centerline
	vtkSmartPointer<vtkConnectivityFilter> connectFilter = vtkSmartPointer<vtkConnectivityFilter>::New();
	connectFilter->SetInputData(m_clippedCenterline);
	connectFilter->SetExtractionModeToAllRegions();
	connectFilter->SetColorRegions(true);
	connectFilter->Update();

	// create implict function with spheres along clipped centerline
	vtkSmartPointer<vtkvmtkPolyBallLine> tubeFunction = vtkSmartPointer<vtkvmtkPolyBallLine>::New();
	tubeFunction->SetInput(m_clippedCenterline);
	tubeFunction->SetPolyBallRadiusArrayName("Radius");

	// loop through all independent patched centerlines
	m_clippedCenterlinesCount = 0;

	vtkNew<vtkImplicitBoolean> endSpheresFunction;
	endSpheresFunction->SetOperationTypeToUnion();

	for (int i = 0; i < connectFilter->GetNumberOfExtractedRegions(); i++)
	{
		// threshold to get independent lines
		vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
		threshold->ThresholdBetween(i, i);
		threshold->SetInputData(connectFilter->GetOutput());
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
		threshold->Update();

		if (threshold->GetOutput()->GetNumberOfPoints() == 0)
			continue;

		// count the number of clipped centerlines
		m_clippedCenterlinesCount++;

		// get the end points
		double* center0 = threshold->GetOutput()->GetPoint(0);
		double tangent0[3];

		tangent0[0] = threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(0,0);
		tangent0[1] = threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(0,1);
		tangent0[2] = threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(0,2);
		double radius0 = threshold->GetOutput()->GetPointData()->GetArray("Radius")->GetComponent(0,0);

		vtkNew<vtkSphere> sphere0;
		sphere0->SetCenter(center0[0], center0[1], center0[2]);
		sphere0->SetRadius(radius0*1.5);

		vtkNew<vtkPlane> plane0;
		plane0->SetOrigin(center0[0], center0[1], center0[2]);
		plane0->SetNormal(1.0*tangent0[0], 1.0*tangent0[1], 1.0*tangent0[2]);

		vtkNew<vtkImplicitBoolean> compositeFunction0;
		compositeFunction0->AddFunction(sphere0);
		compositeFunction0->AddFunction(plane0);
		compositeFunction0->SetOperationTypeToIntersection();

		double* center1 = threshold->GetOutput()->GetPoint(threshold->GetOutput()->GetNumberOfPoints()-1);
		double tangent1[3];
		// reverse tangent direction
		tangent1[0] = -1.0*threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(threshold->GetOutput()->GetNumberOfPoints() - 1, 0);
		tangent1[1] = -1.0*threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(threshold->GetOutput()->GetNumberOfPoints() - 1,1);
		tangent1[2] = -1.0*threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(threshold->GetOutput()->GetNumberOfPoints() - 1,2);
		double radius1 = threshold->GetOutput()->GetPointData()->GetArray("Radius")->GetComponent(0,0);

		vtkNew<vtkSphere> sphere1;
		sphere1->SetCenter(center1[0], center1[1], center1[2]);
		sphere1->SetRadius(radius1*1.5);

		vtkNew<vtkPlane> plane1;
		plane1->SetOrigin(center1[0], center1[1], center1[2]);
		plane1->SetNormal(1.0*tangent1[0], 1.0*tangent1[1], 1.0*tangent1[2]);

		vtkNew<vtkImplicitBoolean> compositeFunction1;
		compositeFunction1->AddFunction(sphere1);
		compositeFunction1->AddFunction(plane1);
		compositeFunction1->SetOperationTypeToIntersection();

		// add to overall composite function
		endSpheresFunction->AddFunction(compositeFunction0);
		endSpheresFunction->AddFunction(compositeFunction1);
	}

	vtkNew<vtkImplicitBoolean> compositeFunction;
	compositeFunction->AddFunction(tubeFunction);
	compositeFunction->AddFunction(endSpheresFunction);
	compositeFunction->SetOperationTypeToDifference();

	// create mask array with spheres and centerline tubes
	vtkSmartPointer<vtkIntArray> maskArray = vtkSmartPointer<vtkIntArray>::New();
	maskArray->SetNumberOfComponents(1);
	maskArray->SetNumberOfTuples(m_voronoiDiagram->GetNumberOfPoints());
	maskArray->SetName("Mask");
	maskArray->FillComponent(0, 0);

	if (m_clippedVoronoiDiagram == NULL)
	{
		m_clippedVoronoiDiagram = vtkSmartPointer<vtkPolyData>::New();
	}
	m_clippedVoronoiDiagram->DeepCopy(m_voronoiDiagram);

	std::cout << "evaluating..." << std::endl;
	compositeFunction->EvaluateFunction(m_clippedVoronoiDiagram->GetPoints()->GetData(),maskArray);

	m_clippedVoronoiDiagram->GetPointData()->AddArray(maskArray);
	m_clippedVoronoiDiagram->GetPointData()->SetActiveScalars("Mask");

	vtkSmartPointer<vtkClipPolyData> clipperV = vtkSmartPointer<vtkClipPolyData>::New();
	clipperV->SetValue(0);
	clipperV->SetInsideOut(true);
	clipperV->GenerateClippedOutputOn();
	clipperV->SetInputData(m_clippedVoronoiDiagram);
	clipperV->Update();

	vtkSmartPointer<vtkCleanPolyData> cleanerV = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanerV->SetInputData(clipperV->GetOutput());
	cleanerV->Update();

	//vtkSmartPointer<vtkGeometryFilter> geomFilter = vtkSmartPointer<vtkGeometryFilter>::New();
	//geomFilter->SetInputData(clipperV->GetOutput());
	//geomFilter->Update();

	m_clippedVoronoiDiagram->DeepCopy(cleanerV->GetOutput());

	std::cout << "save complete" << std::endl;
}

void MouseInteractorStylePickCenterline::InterpoldateVoronoiDiagram()
{
	//// quick load voronoi diagram for debug
	//vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	////reader->SetFileName("Z:\\data\\intracranial\\data_ESASIS_followup\\medical\\055\\baseline\\clipped_voronoi.vtp");
	//reader->SetFileName("D:/Projects/Parent-Vessel-Reconstruction/Data/voronoi_clipped.vtp");
	//reader->Update();
	//m_clippedVoronoiDiagram = vtkSmartPointer<vtkPolyData>::New();
	//m_clippedVoronoiDiagram->DeepCopy(reader->GetOutput());

	////reader->SetFileName("Z:\\data\\intracranial\\data_ESASIS_followup\\medical\\055\\baseline\\centerline_interpolate.vtp");
	//reader->SetFileName("D:/Projects/Parent-Vessel-Reconstruction/Data/centerline_interpolate.vtp");
	//reader->Update();

	//m_centerline->DeepCopy(reader->GetOutput());

	////reader->SetFileName("Z:\\data\\intracranial\\data_ESASIS_followup\\medical\\055\\baseline\\centerline_clipped.vtp");
	//reader->SetFileName("D:/Projects/Parent-Vessel-Reconstruction/Data/centerline_patched.vtp");
	//reader->Update();m_clippedCenterline = vtkSmartPointer<vtkPolyData>::New();;

	//m_clippedCenterline->DeepCopy(reader->GetOutput());

	// interpolate the radius value using clipped voronoi
	vtkNew<vtkPointInterpolator> interpolator;
	interpolator->SetInputData(m_centerline);
	interpolator->SetSourceData(m_clippedVoronoiDiagram);
	interpolator->Update();

	vtkNew<vtkPolyData> newPoints;

	float progress = 0.0;
	const clock_t begin_time = std::clock();

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(m_centerline);
	writer->SetFileName("D:/Projects/Parent-Vessel-Reconstruction/Data/centerline_interpolate.vtp");
	writer->Write();

	writer->SetInputData(m_clippedCenterline);
	writer->SetFileName("D:/Projects/Parent-Vessel-Reconstruction/Data/centerline_patched.vtp");
	writer->Write();

	writer->SetInputData(m_clippedVoronoiDiagram);
	writer->SetFileName("D:/Projects/Parent-Vessel-Reconstruction/Data/voronoi_clipped.vtp");
	writer->Write();

	// loop over independent centerlineids
	for (int i = 0; i < m_clippedCenterline->GetCellData()->GetArray("CenterlineIds")->GetNumberOfTuples(); i++)
	{
		// progress bar
		int barWidth = 30;

		float elapsed_time = float(clock() - begin_time) / CLOCKS_PER_SEC;
		float remain_time = elapsed_time / progress*(1 - progress);

		char et[10];
		char rt[10];
		sprintf(et, "%.2f", elapsed_time);
		sprintf(rt, "%.2f", remain_time);

		std::cout << i << "/" << m_clippedCenterline->GetCellData()->GetArray("CenterlineIds")->GetNumberOfTuples() << "[";
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] " << int(progress * 100.0) << "% "
			<< "Elapsed Time: " << et << "s "
			<< "Remaining Time:" << rt << "s "
			<< "\r";
		std::cout.flush();

		progress = (i + 1)*1.0 / m_clippedCenterline->GetCellData()->GetArray("CenterlineIds")->GetNumberOfTuples();

		// threshold to get independent lines
		vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
		threshold->ThresholdBetween(i, i);
		threshold->SetInputData(m_clippedCenterline);
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "CenterlineIds");
		threshold->Update();

		if (threshold->GetOutput()->GetNumberOfPoints() == 0)
			continue;

		vtkNew<vtkGeometryFilter> geomFilter;
		geomFilter->SetInputData(threshold->GetOutput());
		geomFilter->Update();

		// prepare interpolate input
		vtkPolyData* interpolationEnds = this->ExtractCylindricInterpolationVoronoiDiagram(geomFilter->GetOutput());

		// perform connectivity to get proximal and distal centerlines
		vtkNew<vtkConnectivityFilter> connectFilter;
		connectFilter->SetInputData(geomFilter->GetOutput());
		connectFilter->SetExtractionModeToAllRegions();
		connectFilter->SetColorRegions(true);
		connectFilter->Update();

		if (connectFilter->GetNumberOfExtractedRegions() != 2)
			continue;

		vtkNew<vtkThreshold> thresholdFilter;
		thresholdFilter->SetInputData(connectFilter->GetOutput());
		thresholdFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");

		// get start and end points
		thresholdFilter->ThresholdBetween(0, 0);
		thresholdFilter->Update();
		double startPoint[3];
		startPoint[0] = thresholdFilter->GetOutput()->GetPoint(thresholdFilter->GetOutput()->GetNumberOfPoints()-1)[0];
		startPoint[1] = thresholdFilter->GetOutput()->GetPoint(thresholdFilter->GetOutput()->GetNumberOfPoints() - 1)[1];
		startPoint[2] = thresholdFilter->GetOutput()->GetPoint(thresholdFilter->GetOutput()->GetNumberOfPoints() - 1)[2];

		thresholdFilter->ThresholdBetween(1, 1);
		thresholdFilter->Update();
		double endPoint[3];
		 endPoint[0] = thresholdFilter->GetOutput()->GetPoint(0)[0];
		 endPoint[1] = thresholdFilter->GetOutput()->GetPoint(0)[1];
		 endPoint[2] = thresholdFilter->GetOutput()->GetPoint(0)[2];

		// get start and end point ids
		vtkNew<vtkConnectivityFilter> connectFilter2;
		connectFilter2->SetInputData(interpolator->GetOutput());
		connectFilter2->SetExtractionModeToClosestPointRegion();
		connectFilter2->SetClosestPoint(endPoint);
		connectFilter2->Update();

		vtkNew<vtkPointLocator> centerlineLoactor;
		centerlineLoactor->SetDataSet(connectFilter2->GetOutput());
		centerlineLoactor->BuildLocator();
		int startId = centerlineLoactor->FindClosestPoint(startPoint);
		int endId = centerlineLoactor->FindClosestPoint(endPoint);

		vtkNew<vtkAppendPolyData> appendFilter;
		appendFilter->AddInputData(newPoints);

		// create sphere point clouds
		for (int j = startId; j <= endId; j++)
		{
			// point cloud density update
			int numberOfPoints = connectFilter2->GetOutput()->GetPointData()->GetArray("MaximumInscribedSphereRadius")->GetTuple1(j)*4. / 3.*vtkMath::Pi()*
				pow(connectFilter2->GetOutput()->GetPointData()->GetArray("MaximumInscribedSphereRadius")->GetTuple1(j),3);

			vtkNew<vtkPointSource> pointSource;
			pointSource->SetCenter(connectFilter2->GetOutput()->GetPoint(j));
			pointSource->SetRadius(connectFilter2->GetOutput()->GetPointData()->GetArray("MaximumInscribedSphereRadius")->GetTuple1(j));
			pointSource->SetNumberOfPoints(numberOfPoints); //may need to change to adaptive point cloud density
			pointSource->Update();

			appendFilter->AddInputData(pointSource->GetOutput());
			appendFilter->Update();

			newPoints->DeepCopy(appendFilter->GetOutput());
		}
	}
	std::cout << std::endl;

	// interpolate the radius value to new point cloud
	interpolator->SetInputData(newPoints);
	interpolator->SetSourceData(m_clippedVoronoiDiagram);
	interpolator->Update();

	// merge the new points with the clipped voronoi diagram
	if (m_outputVoronoiDiagram != NULL)
		m_outputVoronoiDiagram->Delete();
	m_outputVoronoiDiagram = vtkSmartPointer<vtkPolyData>::New();

	vtkNew<vtkAppendPolyData> outputAppend;
	outputAppend->AddInputData(m_clippedVoronoiDiagram);
	outputAppend->AddInputDataObject(interpolator->GetOutput());
	outputAppend->Update();

	m_outputVoronoiDiagram->DeepCopy(outputAppend->GetOutput());
}

vtkPolyData* MouseInteractorStylePickCenterline::ExtractCylindricInterpolationVoronoiDiagram(vtkPolyData* clippedCenterline)
{
	std::cout << "Extracting cylindric interpolation Voronoi diagram..." << std::endl;

	// perform connectivity to get proximal and distal centerlines
	vtkNew<vtkConnectivityFilter> connectFilter;
	connectFilter->SetInputData(clippedCenterline);
	connectFilter->SetExtractionModeToAllRegions();
	connectFilter->SetColorRegions(true);
	connectFilter->Update();

	if (connectFilter->GetNumberOfExtractedRegions() != 2)
		return NULL;

	vtkNew<vtkThreshold> thresholdFilter;
	thresholdFilter->SetInputData(connectFilter->GetOutput());
	thresholdFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
	
	// proximal
	thresholdFilter->ThresholdBetween(0, 0);
	thresholdFilter->Update();

	double* tangent1 = thresholdFilter->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetTuple3(thresholdFilter->GetOutput()->GetNumberOfPoints() - 1);
	double tangent1_norm = vtkMath::Norm(tangent1);
	tangent1[0] = tangent1[0] / tangent1_norm;
	tangent1[1] = tangent1[1] / tangent1_norm;
	tangent1[2] = tangent1[2] / tangent1_norm;

	double* point1 = thresholdFilter->GetOutput()->GetPoint(thresholdFilter->GetOutput()->GetNumberOfPoints() - 1);

	// create cylinder function 
	vtkNew<vtkCylinder> cylinderFunction1;
	cylinderFunction1->SetCenter(point1);
	cylinderFunction1->SetRadius(thresholdFilter->GetOutput()->GetPointData()->GetArray("Radius")->GetTuple1(thresholdFilter->GetOutput()->GetNumberOfPoints() - 1)*1.5);
	cylinderFunction1->SetAxis(tangent1);

	// clip planes
	vtkNew<vtkPlane> planeFunction11;
	double plane11origin[3];
	plane11origin[0] = point1[0] - tangent1[0]*1.5;
	plane11origin[1] = point1[1] - tangent1[1] * 1.5;
	plane11origin[2] = point1[2] - tangent1[2] * 1.5;

	planeFunction11->SetOrigin(plane11origin);
	planeFunction11->SetNormal(tangent1[0] * -1, tangent1[1] * -1, tangent1[2] * -1);

	vtkNew<vtkPlane> planeFunction12;
	double plane12origin[3];
	plane12origin[0] = point1[0] + tangent1[0] * 1.5;
	plane12origin[1] = point1[1] + tangent1[1] * 1.5;
	plane12origin[2] = point1[2] + tangent1[2] * 1.5;

	planeFunction12->SetOrigin(plane12origin);
	planeFunction12->SetNormal(tangent1[0]*1, tangent1[1] * 1, tangent1[2] * 1);
	
	// distal 
	thresholdFilter->ThresholdBetween(1, 1);
	thresholdFilter->Update();

	double* tangent2 = thresholdFilter->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetTuple3(0);
	double tangent2_norm = vtkMath::Norm(tangent2);
	tangent2[0] = tangent2[0] / tangent2_norm;
	tangent2[1] = tangent2[1] / tangent2_norm;
	tangent2[2] = tangent2[2] / tangent2_norm;

	double* point2 = thresholdFilter->GetOutput()->GetPoint(0);

	// create cylinder function 
	vtkNew<vtkCylinder> cylinderFunction2;
	cylinderFunction2->SetCenter(point2);
	cylinderFunction2->SetRadius(thresholdFilter->GetOutput()->GetPointData()->GetArray("Radius")->GetTuple1(0) *1.5);
	cylinderFunction2->SetAxis(tangent2);

	// clip planes
	vtkNew<vtkPlane> planeFunction21;
	double plane21origin[3];
	plane21origin[0] = point2[0] - tangent2[0] * 0.5;
	plane21origin[1] = point2[1] - tangent2[1] * 0.5;
	plane21origin[2] = point2[2] - tangent2[2] * 0.5;

	planeFunction21->SetOrigin(plane21origin);
	planeFunction21->SetNormal(tangent2[0] * -1, tangent2[1] * -1, tangent2[2] * -1);

	vtkNew<vtkPlane> planeFunction22;
	double plane22origin[3];
	plane22origin[0] = point2[0] + tangent2[0] * 0.5;
	plane22origin[1] = point2[1] + tangent2[1] * 0.5;
	plane22origin[2] = point2[2] + tangent2[2] * 0.5;

	planeFunction22->SetOrigin(plane22origin);
	planeFunction22->SetNormal(tangent2[0] * 1, tangent2[1] * 1, tangent2[2] * 1);

	// evaluate the points
	// convert clipped voronoi points to vtk data array
	vtkNew<vtkDoubleArray> mask;
	mask->SetNumberOfComponents(1);
	mask->SetNumberOfTuples(m_clippedVoronoiDiagram->GetNumberOfPoints());
	mask->SetName("Mask");
	mask->FillComponent(0, 0);

	vtkNew<vtkImplicitBoolean> composeFunction1;
	composeFunction1->AddFunction(planeFunction11);
	composeFunction1->AddFunction(planeFunction12);
	composeFunction1->AddFunction(cylinderFunction1);
	composeFunction1->SetOperationTypeToIntersection();
	composeFunction1->EvaluateFunction(m_clippedVoronoiDiagram->GetPoints()->GetData(), mask);

	vtkNew<vtkPolyData> interpolationEnd1;
	interpolationEnd1->DeepCopy(m_clippedVoronoiDiagram);
	interpolationEnd1->GetPointData()->AddArray(mask);
	interpolationEnd1->GetPointData()->SetActiveScalars("Mask");

	vtkNew<vtkClipPolyData> clipper1;
	clipper1->SetValue(1e-6);
	clipper1->SetInsideOut(true);
	clipper1->GenerateClippedOutputOn();
	clipper1->SetInputData(interpolationEnd1);
	clipper1->Update();

	vtkSmartPointer<vtkCleanPolyData> cleaner1 = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner1->SetInputData(clipper1->GetOutput());
	cleaner1->Update();

	double vol1 =
		(cleaner1->GetOutput()->GetBounds()[1] - cleaner1->GetOutput()->GetBounds()[0])*
		(cleaner1->GetOutput()->GetBounds()[3] - cleaner1->GetOutput()->GetBounds()[2])*
		(cleaner1->GetOutput()->GetBounds()[5] - cleaner1->GetOutput()->GetBounds()[4]);
	double density1 = cleaner1->GetOutput()->GetNumberOfPoints() / vol1;

	// distal end
	vtkNew<vtkImplicitBoolean> composeFunction2;
	composeFunction2->AddFunction(planeFunction21);
	composeFunction2->AddFunction(planeFunction22);
	composeFunction2->AddFunction(cylinderFunction2);
	composeFunction2->SetOperationTypeToIntersection();
	composeFunction2->EvaluateFunction(m_clippedVoronoiDiagram->GetPoints()->GetData(), mask);

	vtkNew<vtkPolyData> interpolationEnd2;
	interpolationEnd2->DeepCopy(m_clippedVoronoiDiagram);
	interpolationEnd2->GetPointData()->AddArray(mask);
	interpolationEnd2->GetPointData()->SetActiveScalars("Mask");

	vtkNew<vtkClipPolyData> clipper2;
	clipper2->SetValue(1e-6);
	clipper2->SetInsideOut(true);
	clipper2->GenerateClippedOutputOn();
	clipper2->SetInputData(interpolationEnd2);
	clipper2->Update();

	vtkSmartPointer<vtkCleanPolyData> cleaner2 = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner2->SetInputData(clipper2->GetOutput());
	cleaner2->Update();

	double vol2 =
		(cleaner2->GetOutput()->GetBounds()[1] - cleaner2->GetOutput()->GetBounds()[0])*
		(cleaner2->GetOutput()->GetBounds()[3] - cleaner2->GetOutput()->GetBounds()[2])*
		(cleaner2->GetOutput()->GetBounds()[5] - cleaner2->GetOutput()->GetBounds()[4]);
	double density2 = cleaner2->GetOutput()->GetNumberOfPoints() / vol2;

	// square to sphere
	double density_correction_factor = 1./(vtkMath::Pi() / 6.);

	m_pointCloudDensity = (density1 + density2) / 2*density_correction_factor;

	// combine the masks
	vtkNew<vtkImplicitBoolean> composeFunction3;
	composeFunction3->AddFunction(composeFunction1);
	composeFunction3->AddFunction(composeFunction2);
	composeFunction3->SetOperationTypeToUnion();
	composeFunction3->EvaluateFunction(m_clippedVoronoiDiagram->GetPoints()->GetData(), mask);
	
	vtkNew<vtkPolyData> interpolationEnds;
	interpolationEnds->DeepCopy(m_clippedVoronoiDiagram);
	interpolationEnds->GetPointData()->AddArray(mask);

	return interpolationEnds;
}

vtkStandardNewMacro(MouseInteractorStylePickCenterline);
