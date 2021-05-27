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

#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"

#include "vtkvmtkPolyBallLine.h"

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
		this->ClipVoronoiDiagram();
		//this->InterpoldateVoronoiDiagram();
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

	std::vector<std::pair<vtkSphere*,double*>> endPointSpheres;

	// loop through all independent patched centerlines
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

		// get the end points
		double* center0 = threshold->GetOutput()->GetPoint(0);
		double* tangent0 = new double [3];

		tangent0[0] = threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(0,0);
		tangent0[1] = threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(0,1);
		tangent0[2] = threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(0,2);
		double radius0 = threshold->GetOutput()->GetPointData()->GetArray("Radius")->GetComponent(0,0);

		vtkSphere* sphere0 = vtkSphere::New();
		sphere0->SetCenter(center0[0], center0[1], center0[2]);
		sphere0->SetRadius(radius0*1.5);

		std::pair<vtkSphere*, double*> endPoint0;
		endPoint0.first = sphere0;
		endPoint0.second = tangent0;
		endPointSpheres.push_back(endPoint0);

		double* center1 = threshold->GetOutput()->GetPoint(threshold->GetOutput()->GetNumberOfPoints()-1);
		double* tangent1 = new double [3];
		// reverse tangent direction
		tangent1[0] = -1.0*threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(threshold->GetOutput()->GetNumberOfPoints() - 1, 0);
		tangent1[1] = -1.0*threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(threshold->GetOutput()->GetNumberOfPoints() - 1,1);
		tangent1[2] = -1.0*threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetComponent(threshold->GetOutput()->GetNumberOfPoints() - 1,2);
		double radius1 = threshold->GetOutput()->GetPointData()->GetArray("Radius")->GetComponent(0,0);

		vtkSphere* sphere1 = vtkSphere::New();

		sphere1->SetCenter(center1[0], center1[1], center1[2]);
		sphere1->SetRadius(radius1*1.5);

		std::pair<vtkSphere*, double*> endPoint1;
		endPoint1.first = sphere1;
		endPoint1.second = tangent1;
		endPointSpheres.push_back(endPoint1);
	}

	// create mask array with spheres and centerline tubes
	vtkSmartPointer<vtkIntArray> maskArray = vtkSmartPointer<vtkIntArray>::New();
	maskArray->SetNumberOfComponents(1);
	maskArray->SetNumberOfTuples(m_voronoiDiagram->GetNumberOfPoints());
	maskArray->SetName("Mask");
	maskArray->FillComponent(0, 0);

	float progress = 0.0;
	const clock_t begin_time = std::clock();

	for (int i = 0; i < m_voronoiDiagram->GetNumberOfPoints(); i++)
	{
		double* point = m_voronoiDiagram->GetPoint(i);

		// loop over end spheres
		bool in_end_point_sphere_and_oppose_tangent = false;

		for (int j = 0; j < endPointSpheres.size(); j++)
		{
			double center2voronoi[3];
			double* c = endPointSpheres.at(j).first->GetCenter();

			center2voronoi[0] = point[0] - c[0];
			center2voronoi[1] = point[1] - c[1];
			center2voronoi[2] = point[2] - c[2];
			
			double directionDot = vtkMath::Dot(center2voronoi, endPointSpheres.at(j).second);

			if (endPointSpheres.at(j).first->EvaluateFunction(point) < 0.0 && directionDot < 0.0)
			{
				in_end_point_sphere_and_oppose_tangent = true;
				break;
			}
		}

		if (in_end_point_sphere_and_oppose_tangent)
		{
			continue;
		}
			
		double tubeValue = tubeFunction->EvaluateFunction(point);
		if (tubeValue <= 0)
		{
			maskArray->SetTuple1(i, 1);
		}

		// progress bar
		int barWidth = 30;

		float elapsed_time = float(clock() - begin_time) / CLOCKS_PER_SEC;
		float remain_time = elapsed_time / progress*(1 - progress);

		char et[10];
		char rt[10];
		sprintf(et, "%.2f", elapsed_time);
		sprintf(rt, "%.2f", remain_time);

		std::cout << i << "/" << m_voronoiDiagram->GetNumberOfPoints() <<"[";
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

		progress = (i+1)*1.0/ m_voronoiDiagram->GetNumberOfPoints();
	}
	std::cout << std::endl;

	if (m_outputVoronoiDiagram == NULL)
	{
		m_outputVoronoiDiagram = vtkSmartPointer<vtkPolyData>::New();
	}
	m_outputVoronoiDiagram->DeepCopy(m_voronoiDiagram);
	m_outputVoronoiDiagram->GetPointData()->AddArray(maskArray);
	m_outputVoronoiDiagram->GetPointData()->SetActiveScalars("Mask");

	vtkSmartPointer<vtkClipPolyData> clipperV = vtkSmartPointer<vtkClipPolyData>::New();
	clipperV->SetValue(1e-6);
	clipperV->SetInsideOut(false);
	clipperV->GenerateClippedOutputOn();
	clipperV->SetInputData(m_outputVoronoiDiagram);
	clipperV->Update();

	vtkSmartPointer<vtkCleanPolyData> cleanerV = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanerV->SetInputData(clipperV->GetOutput());
	cleanerV->Update();

	vtkSmartPointer<vtkGeometryFilter> geomFilter = vtkSmartPointer<vtkGeometryFilter>::New();
	geomFilter->SetInputData(cleanerV->GetOutput());
	geomFilter->Update();

	if (m_clippedVoronoiDiagram == NULL)
	{
		m_clippedVoronoiDiagram = vtkSmartPointer<vtkPolyData>::New();
	}

	m_clippedVoronoiDiagram->DeepCopy(clipperV->GetOutput());

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(m_centerline);
	writer->SetFileName("Z:/data/intracranial/data_ESASIS_followup/medical/055/baseline/centerline_interpolate.vtp");
	writer->Write();

	writer->SetInputData(m_clippedCenterline);
	writer->SetFileName("Z:/data/intracranial/data_ESASIS_followup/medical/055/baseline/centerline_patched.vtp");
	writer->Write();

	writer->SetInputData(m_outputVoronoiDiagram);
	writer->SetFileName("Z:/data/intracranial/data_ESASIS_followup/medical/055/baseline/voronoi_masked.vtp");
	writer->Write();

	writer->SetInputData(m_clippedVoronoiDiagram);
	writer->SetFileName("Z:/data/intracranial/data_ESASIS_followup/medical/055/baseline/voronoi_clipped.vtp");
	writer->Write();

	std::cout << "save complete" << std::endl;

	//recycle the spheres
	while (endPointSpheres.size() > 0)
	{
		endPointSpheres.back().first->Delete();
		delete endPointSpheres.back().second;
		endPointSpheres.pop_back();
	}
}

void MouseInteractorStylePickCenterline::InterpoldateVoronoiDiagram()
{
	// quick load voronoi diagram for debug
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName("Z:\\data\\intracranial\\data_ESASIS_followup\\medical\\055\\baseline\\clipped_voronoi.vtp");
	reader->Update();
	m_clippedVoronoiDiagram->DeepCopy(reader->GetOutput());

	reader->SetFileName("Z:\\data\\intracranial\\data_ESASIS_followup\\medical\\055\\baseline\\normalized_centerline.vtp");
	reader->Update();
	m_centerline->DeepCopy(reader->GetOutput());

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(m_clippedVoronoiDiagram);
	writer->SetFileName("Z:/data/intracranial/data_ESASIS_followup/medical/055/baseline/clipped_voronoi.vtp");
	writer->Write();
	std::cout << "save complete" << std::endl;
}

vtkStandardNewMacro(MouseInteractorStylePickCenterline);
