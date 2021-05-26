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

		//m_centerline->DeepCopy(converter->GetOutput());
		//vtkSmartPointer <vtkXMLPolyDataWriter> writer = vtkSmartPointer <vtkXMLPolyDataWriter>::New();
		//writer->SetInputData(clipper2->GetOutput());
		//writer->SetFileName("Z:/data/intracranial/data_ESASIS_followup/medical/055/baseline/normalized_centerline.vtp");
		//writer->Write();

		//std::cout << "save ok" << std::endl;
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

			//std::cout << "picked point: " <<
			//	m_centerline->GetPoint(m_pickedPointId)[0] << ", " << 
			//	m_centerline->GetPoint(m_pickedPointId)[1] << ", " <<
			//	m_centerline->GetPoint(m_pickedPointId)[2] << std::endl;

			//std::cout << "com: " <<
			//	comFilter->GetCenter()[0] << ", " <<
			//	comFilter->GetCenter()[1] << ", " <<
			//	comFilter->GetCenter()[2] << std::endl;

			double centroid_dist =
				vtkMath::Distance2BetweenPoints(
					m_centerline->GetPoint(m_pickedPointId),
					comFilter->GetCenter()
				);

			//std::cout << "centroid_dist: " << centroid_dist << std::endl;

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
		
		this->ClipVoronoiDiagram();

		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetInputData(m_clippedCenterline);
		writer->SetFileName("Z:/data/intracranial/data_ESASIS_followup/medical/055/baseline/normalized_centerline.vtp");
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

		std::cout << "cetnerlineid:" << i << ", number of points: " << threshold->GetOutput()->GetNumberOfPoints() << std::endl;

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

	//std::cout << "origin: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;

	double point1[3];
	point1[0] = origin[0] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id)[0] / normal_norm;
	point1[1] = origin[1] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id)[1] / normal_norm;
	point1[2] = origin[2] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetNormal")->GetTuple(id)[2] / normal_norm;

	//std::cout << "point1: " << point1[0] << " " << point1[1] << " " << point1[2] << std::endl;

	double point2[3];
	point2[0] = origin[0] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id)[0] / binormal_norm;
	point2[1] = origin[1] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id)[1] / binormal_norm;
	point2[2] = origin[2] + sqrt(2)*size_factor * radius*centerline->GetPointData()->GetArray("FrenetBinormal")->GetTuple(id)[2] / binormal_norm;
	//std::cout << "point2: " << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;

	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
	plane->SetOrigin(origin);
	plane->SetPoint1(point1);
	plane->SetPoint2(point2);
	plane->Update();

	//std::cout << "create plane ok" << std::endl;

	vtkSmartPointer<vtkPolyDataMapper> priximal_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	priximal_mapper->SetInputData(plane->GetOutput());

	actor->SetMapper(priximal_mapper);
	actor->GetProperty()->SetColor(1, 0, 0);
}

void MouseInteractorStylePickCenterline::ClipVoronoiDiagram()
{
	//if (m_voronoiDiagram == NULL)
	//	return;

	vtkSmartPointer<vtkIntArray> maskArray = vtkSmartPointer<vtkIntArray>::New();
	maskArray->SetNumberOfComponents(1);
	maskArray->SetNumberOfTuples(m_voronoiDiagram->GetNumberOfPoints());
	maskArray->FillComponent(0, 0);

	// compute connectivity on clipped centerline
	vtkSmartPointer<vtkConnectivityFilter> connectFilter = vtkSmartPointer<vtkConnectivityFilter>::New();
	connectFilter->SetInputData(m_clippedCenterline);
	connectFilter->SetExtractionModeToAllRegions();
	connectFilter->Update();

	// create implict function with spheres along clipped centerline
	vtkSmartPointer<vtkvmtkPolyBallLine> tubeFunction = vtkSmartPointer<vtkvmtkPolyBallLine>::New();
	tubeFunction->SetInput(m_clippedCenterline);
	tubeFunction->SetPolyBallRadiusArrayName("Radius");

	std::list<vtkSphere*> 

	// loop through all independent patched centerlines
	for (int i = 0; i < connectFilter->GetNumberOfExtractedRegions(); i++)
	{
		// threshold to get independent lines
		vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
		threshold->ThresholdBetween(i, i);
		threshold->SetInputData(connectFilter->GetOutput());
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
		threshold->Update();

		// get the end points
		double* center0 = threshold->GetOutput()->GetPoint(0);
		double* tangent0 = threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetTuple(0);

		double* center1 = threshold->GetOutput()->GetPoint(threshold->GetOutput()->GetNumberOfPoints()-1);
		double* tangent1 = threshold->GetOutput()->GetPointData()->GetArray("FrenetTangent")->GetTuple(threshold->GetOutput()->GetNumberOfPoints() - 1);
		// reverse tangent direction
		tangent1[0] = -1.0*tangent1[0];
		tangent1[1] = -1.0*tangent1[1];
		tangent1[2] = -1.0*tangent1[2];
	}



	//	if (threshold->GetOutput()->GetNumberOfPoints() == 0)
	//		continue;

	//	// compute patch end point parameters

	//	double point0[3];
	//	double point1[3];
	//	double radius0;

	//	if (i == 0)
	//	{
	//		cell->GetPoints()->GetPoint(cell->GetNumberOfPoints() - 1,point0);
	//		cell->GetPoints()->GetPoint(cell->GetNumberOfPoints() - 2, point1);
	//		radius0 = patchCenterline->GetPointData()->GetArray("Radius")
	//			->GetTuple1(cell->GetPointId(cell->GetNumberOfPoints() - 1));
	//	}
	//	else
	//	{
	//		cell->GetPoints()->GetPoint(0, point0);
	//		cell->GetPoints()->GetPoint(1, point1);
	//		radius0 = patchCenterline->GetPointData()->GetArray("Radius")
	//			->GetTuple1(cell->GetPointId(0));
	//	}

	//	double tan[3];
	//	tan[0] = point1[0] - point0[0];
	//	tan[1] = point1[1] - point0[1];
	//	tan[2] = point1[2] - point0[2];
	//	vtkMath::Normalize(tan);

	//	return tan, point0, radius0

	//	// mask with patch
	//}
}

vtkStandardNewMacro(MouseInteractorStylePickCenterline);
