#include <QApplication>
#include <QCommandLineOption>
#include <QCommandLineParser>

#include "ParentVesselReconstruction.h"

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

	ParentVesselReconstruction pvr;
	pvr.SetSourceFilePath(args.at(0));
	pvr.SetOutputFilePath(args.at(1));
	pvr.Run();

	return app.exec();
}