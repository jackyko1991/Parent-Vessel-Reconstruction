#include <QApplication>
#include <QCommandLineOption>
#include <QCommandLineParser>
#include <QDir>

#include "ParentVesselReconstruction.h"

#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "observe_error.h""
#include "vtkSTLReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"

#include <iostream>

int main(int argc, char** argv)
{
	QApplication app(argc, argv);
	QApplication::setApplicationName("Parent Vessel Reconstruction");
	QApplication::setApplicationVersion("1.0");

	QCommandLineParser parser;
	parser.setApplicationDescription("Reconstruct parent vessel from aneurysm vessels");
	parser.addHelpOption();
	parser.addVersionOption();
	parser.addPositionalArgument("source", QCoreApplication::translate("main", "Source file to process."));
	parser.addPositionalArgument("centerline", QCoreApplication::translate("main", "Centerline file to process."));
	parser.addPositionalArgument("output", QCoreApplication::translate("main", "Output file location."));
	
	QCommandLineOption intermediateFilesOption("i", QCoreApplication::translate("main", "Save intermediate files"));
	parser.addOption(intermediateFilesOption);

	parser.process(app);
	const QStringList args = parser.positionalArguments();

	if (args.size() < 2)
	{
		std::cerr << "Source/output file not provided. Use \"-h\" to see usage." << std::endl;
		return 0;
	}

	QFileInfo sourceFile(args[0]);
	QFileInfo centerlineFile(args[1]);

	// check input file existence
	if (!sourceFile.exists())
	{
		std::cout << "Source file not exist!!!" << std::endl;
		return 0;
	}

	// read surface
	vtkSmartPointer<ErrorObserver> errorObserver = vtkSmartPointer<ErrorObserver>::New();

	vtkNew<vtkPolyData> source;
	vtkNew<vtkPolyData> centerline;

	if (sourceFile.suffix() == "vtp" || sourceFile.suffix() == "VTP")
	{
		vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
		reader->SetFileName(sourceFile.absoluteFilePath().toStdString().c_str());
		reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
		reader->Update();
		source->DeepCopy(reader->GetOutput());
	}
	else if (sourceFile.suffix() == "stl" || sourceFile.suffix() == "STL")
	{
		vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(sourceFile.absoluteFilePath().toStdString().c_str());
		reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
		reader->Update();
		source->DeepCopy(reader->GetOutput());

	}
	else if (sourceFile.suffix() == "vtk" || sourceFile.suffix() == "VTK")
	{
		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
		reader->SetFileName(sourceFile.absoluteFilePath().toStdString().c_str());
		reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
		source->DeepCopy(reader->GetOutput());
	}
	else
	{
		std::cerr << "Invalid input data type, only accept STL, VTP or VTK files" << std::endl;
		return 0;
	}

	std::cout << "Read source file complete" << std::endl;

	// read centerline
	if (centerlineFile.suffix() == "vtp" || centerlineFile.suffix() == "VTP")
	{
		vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
		reader->SetFileName(centerlineFile.absoluteFilePath().toStdString().c_str());
		reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
		reader->Update();
		centerline->DeepCopy(reader->GetOutput());
	}
	else if (centerlineFile.suffix() == "vtk" || centerlineFile.suffix() == "VTK")
	{
		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
		reader->SetFileName(centerlineFile.absoluteFilePath().toStdString().c_str());
		reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
		centerline->DeepCopy(reader->GetOutput());
	}
	else
	{
		std::cerr << "Invalid centerline data type, only accept VTP or VTK files" << std::endl;
		return 0;
	}

	std::cout << "Read centerline complete" << std::endl;

	ParentVesselReconstruction pvr;
	pvr.SetSource(source);
	pvr.SetCenterline(centerline);
	pvr.Run();

	// the 3rd argument is save director
	QDir outputDir;
	if (args.size() > 2)
		outputDir.setPath(args[2]);
	else
		return 0;

	// check if Directory exists
	if (!outputDir.exists())
	{
		// create output directory
		outputDir.mkpath(outputDir.absolutePath());
	}

	std::cout << "Writing output files to " << outputDir.absolutePath().toStdString() << std::endl;

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	QString outputName;
	writer->SetInputData(pvr.GetVoronoiDiagram());
	outputName = outputDir.absolutePath() + "/voronoi_recon.vtp";
	writer->SetFileName(outputName.toStdString().c_str());
	writer->Write();

	writer->SetInputData(pvr.GetClippedVoronoiDiagram());
	outputName = outputDir.absolutePath() + "/voronoi_clipped.vtp";
	writer->SetFileName(outputName.toStdString().c_str());
	writer->Write();

	writer->SetInputData(pvr.GetCenterline());
	outputName = outputDir.absolutePath() + "/centerline_interpolate.vtp";
	writer->SetFileName(outputName.toStdString().c_str());
	writer->Write();

	writer->SetInputData(pvr.GetClippedCenterline());
	outputName = outputDir.absolutePath() + "/centerline_patched.vtp";

	writer->SetFileName(outputName.toStdString().c_str());
	writer->Write();
	//
	//writer->SetInputData(pvr.GetSource());
	//writer->SetFileName("D:/Projects/Parent-Vessel-Reconstruction/Data/surface_reconstructed.vtp");
	//writer->Write();

	std::cout << "save complete" << std::endl;
}