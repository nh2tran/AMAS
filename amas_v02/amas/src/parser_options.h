// =====================================================================
// Class "Options" and functions "parse..."
// =====================================================================






#ifndef NHTRAN_PARSER_OPTIONS_H_
#define NHTRAN_PARSER_OPTIONS_H_






#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>






using namespace seqan;






// =====================================================================
// Class "Options"
// =====================================================================

struct Options {
	
	CharString genomeFile;
	CharString readsFile;
	CharString readsFile2;

	//unsigned F; // frequency threshold of adaptive seeds
	
	unsigned errors; // maximum number of errors allowed (mismatches & gaps)
	
	//unsigned msl; // maximum seed length
	
	unsigned fse; // require full-sensitivity for matches up to "fse" errors
	
	bool best; // Best-mapping mode
	
	unsigned k; // k-mapping mode
	
	bool raw; // output the alignments in the raw format for pair-end mapping with "masai_output_pe"
	
	unsigned threads; // number of threads allocated
	
	// Debug Options - same as those in MASAI
	// to record the seeding time or the mapping time only
    bool noVerify; // seed partition only, do not perform seed extension
    bool noDump; // mapping only, do not write alignment results

	Options():
		errors(5),
		fse(2),
		best(false),
		k(1000000),
		raw(false),
		threads(1),
		noVerify(false),
		noDump(false)
	{}
};











// =====================================================================
// Function parseHasher()
// =====================================================================

ArgumentParser::ParseResult 
parseHasher(Options & options, int argc, char const ** argv) 
{
	
	ArgumentParser parser; // defined in <seqan/arg_parse.h>
	
    setAppName(parser, "amas_indexer");
    setShortDescription(parser, "Build the adaptive-seed tree from the reference genome");
    setVersion(parser, "0.2");
    setDate(parser, "February 2014");

    addUsageLine(parser, "<\\fIGENOME FILE\\fP>");
	//addUsageLine(parser, "<\\fIGENOME FILE\\fP> [\\fIOPTIONS\\fP]");

    addDescription(parser, "Please provide the genome file in FASTA format with extension \".fa\" or \".fasta\".");
    //addDescription(parser, "The FM index of the genome must be built beforehand (using amas_indexer).");
    addDescription(parser, "The output tree will be stored in folder \"tree\".");
    addDescription(parser, "The full name of the genome file will be used as the prefix of the tree file's name.");
    addDescription(parser, "This app has no parameters for configuration by users.");

    
    
    
    
    
    // ARGUMENTS
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta fa");
    
    
    
    
    

    // parse command line
    ArgumentParser::ParseResult parseResult = parse(parser, argc, argv);
    
    if (parseResult != ArgumentParser::PARSE_OK) {
		return parseResult;
	}
	
	
	
	
	
	
	// extract arguments & options from "parser"
	getArgumentValue(options.genomeFile, parser, 0);
	
	
	
	
	
	
    return ArgumentParser::PARSE_OK;
	
}











// =====================================================================
// Function parseMapper()
// =====================================================================
ArgumentParser::ParseResult 
parseMapper(Options & options, int argc, char const ** argv) 
{
	
	ArgumentParser parser; /* defined in <seqan/arg_parse.h> */
	
    setAppName(parser, "amas_mapper");
    setShortDescription(parser, "An NGS all-mapping tool using adaptive seeds");
    setVersion(parser, "0.2");
    setDate(parser, "February 2014");

    addUsageLine(parser, "<\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP> [\\fIOPTIONS\\fP]");

    addDescription(parser, "Please provide the genome file in FASTA format with extension \".fa\" or \".fasta\".");
    addDescription(parser, "The reads file must be in FASTA or FASTQ format with extension \".fa\", \".fasta\", \".fq\", or \".fastq\".");
    addDescription(parser, "The adaptive-seed tree of the genome must be built beforehand (using amas_indexer).");
    addDescription(parser, "The output file is named as following \"reads_file.amas.sam\".");
    
    
    
    
    
    
    // ARGUMENTS
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta fa");
    setValidValues(parser, 1, "fasta fa fastq fq");
    




	
	// OPTIONS
    addSection(parser, "OPTIONS");
    
    addOption(parser, ArgParseOption("e", "errors", "Maximum number of errors allowed (mismatches & gaps).", ArgParseOption::INTEGER));
    setMinValue(parser, "errors", "0");
    setMaxValue(parser, "errors", "10");
    setDefaultValue(parser, "errors", 5);

    addOption(parser, ArgParseOption("fse", "full-sen-error", "Require full sensitivity for matches with up to \"fse\" errors.", ArgParseOption::INTEGER));
    setMinValue(parser, "full-sen-error", "0");
    setMaxValue(parser, "full-sen-error", "5");
    setDefaultValue(parser, "full-sen-error", 2);

    addOption(parser, ArgParseOption("best", "best", "Best-mapping mode: find the best alignment for each read."));

    addOption(parser, ArgParseOption("k", "k-mapping", "k-mapping mode: find up to k alignments for each reads.", ArgParseOption::INTEGER));
    setMinValue(parser, "k-mapping", "1");
    setMaxValue(parser, "k-mapping", "1000000");
    setDefaultValue(parser, "k-mapping", 1000000);

    addOption(parser, ArgParseOption("raw", "raw", "Output the alignments in the raw format for pair-end matching with \"masai_output_pe\"."));

    addOption(parser, ArgParseOption("p", "threads", "Number of threads.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "16");
    setDefaultValue(parser, "threads", 1);






	// Debug Options
    //addSection(parser, "Debug Options");
    //addOption(parser, ArgParseOption("nv", "no-verify", "Do not perform seed extension."));
    //addOption(parser, ArgParseOption("nd", "no-dump", "Do not write alignment results."));






    // parse command line
    ArgumentParser::ParseResult parseResult = parse(parser, argc, argv);
    
    if (parseResult != ArgumentParser::PARSE_OK) {
		return parseResult;
	}
	
	
	
	
	
	
	// extract arguments & options from "parser"
	getArgumentValue(options.genomeFile, parser, 0);
	getArgumentValue(options.readsFile, parser, 1);
	getOptionValue(options.errors, parser, "errors");
	getOptionValue(options.fse, parser, "full-sen-error");
	options.best = isSet(parser, "best");
	getOptionValue(options.k, parser, "k-mapping");
    options.raw = isSet(parser, "raw");
	getOptionValue(options.threads, parser, "threads");
	
    //options.noVerify = isSet(parser, "no-verify");
    //options.noDump = isSet(parser, "no-dump");
    
    
    
    
	
	
    return ArgumentParser::PARSE_OK;
	
}

#endif
