#!/usr/bin/perl
use 5.010;
use strict;
use warnings;
use Term::Prompt;

## This program is to check the existing modules in your system. You can add more module in @modules array 

my $script = 'EBA.pl';

die "Provide script file name on the command line\n"
    unless defined $script;

	
my @modules = qw(
	File::Path
	Math::Round
	DateTime::Locale
	Term::Prompt
	Term::ANSIColor
	List::Compare
    	List::Compare::Base::_Engine
   	List::Compare::Base::_Auxiliary
);

print "Checking mandatory modules for EBA\n";
checkModules(@modules); ## Provide the name of all modules
	
until ( do $script ) {
    my $expire = $@;
    if ( my ($file) = $expire =~ /^Can't locate (.+?) in/ ) {
        my $module = $file;
        $module =~ s/\.(\w+)$//;
        $module = join('::', split '/', $module);
        print "Attempting to install '$module' via cpan\n";
        system(cpan => $module);
        last unless prompt(y => 'Try Again?', '', 'n');
    }
    else {
        die $expire;
    }
}


sub checkModules {
	my @reqMod=@_;
	for(@reqMod) {
    	eval "use $_";
    		if ($@) {
        		warn "Not found -> $_ module\t Need to install $_ \n" if $@;
    		} else {
        		say "Found : $_";
    		}
	}
}



__END__

## List of all modules used in EBA tool

'drawBetaGraph'                   => 'undef',
'drawBreakpointGraph'             => 'undef',
'BreaksAmongstSpecies'            => 'undef',
'BreaksFinder'                    => 'undef',
'BreaksMatrix'                    => 'undef',
'BreaksScoring'                   => 'undef',
'BreaksScoring2'                  => 'undef',
'CalculateBeta'                   => 'undef',
'CheckData'                       => 'undef',
'ClassifyBreakpoints'             => 'undef',
'ConCatAll'                       => 'undef',
'ConCatFile'                      => 'undef',
'CreateFinal'                     => 'undef',
'EnterScore'                      => 'undef',
'FindReuse'                       => 'undef',
'FindReuseMerge'                  => 'undef',
'MergeResolution'                 => 'undef',
'ModifyBeta'                      => 'undef',
'ModifyClassification'            => 'undef',
'Safai'                           => 'undef',
'StoreSpecies'                    => 'undef',
'Visualize'                       => 'undef',
'drawBreakpointChrGraphFinal'     => 'undef',
'drawBrkChrLineGraphFinal'        => 'undef',
'drawCumulatedBar'                => 'undef',
'drawCumulatedStackedBar'         => 'undef',
'drawCumulatedStackedBarAll'      => 'undef',
'drawCumulatedStackedBarMerged'   => 'undef',
'drawPieChart'                    => 'undef',
'drawPieChartBreaksFinal'         => 'undef',
'drawPieChartBreaksFinalMerged'   => 'undef',
'drawPieChartFinal'               => 'undef',
'drawChrBreakpointGraph'          => 'undef',
'GanchoApproach'                  => 'undef',
'List::Compare'                   => '0.37',
'GD::Polygon'                     => 'undef',
'GD::Graph::Data'                 => '1.22',
'GD::Image'                       => '2.38',
'GD::Graph::Error'                => '1.8',
'GD::Text'                        => '0.86',
'GD::Graph::bars'                 => '1.26',
'GD::Graph::axestype'             => '1.45',
'GD'                              => '2.46',
'GD::Text::Align'                 => '1.18',
'GD::Graph::colour'               => '1.10',
'GD::Graph::utils'                => '1.7',
'Math::Round'                     => '0.05',
'List::Compare::Base::_Engine'    => '0.37',
'List::Compare::Base::_Auxiliary' => '0.37',
'GD::Graph::lines'                => '1.15',
'GD::Graph::hbars'                => '1.3',
'GD::Graph'                       => '1.44',
'GD::Graph::pie'                  => '1.21',


##To Install Modules

#!/usr/bin/perl

## Script to install perl module .. 
use strict;
use warnings;
use CPAN;

CPAN::Shell->install(

"Test::PDF",
"ExtUtils::PkgConfig",
"ExtUtils::Depends",
"Cairo",
"Test::LongString",
"Text::Flow",
"Graphics::Primitive::Driver::Cairo",
"Math::Gradient",
"Set::Infinite",
"DateTime::Set",
"Algorithm::Diff",
"Text::Diff",
"Test::Differences",
"Color::Scheme",
"Class::Data::Inheritable",
"Class::Accessor::Fast",
"Module::Pluggable",
"Color::Library",
"Test::Number::Delta",
"MooseX::Aliases",
"Graphics::Color",
"Forest",
"Graphics::Primitive",
"MooseX::AttributeHelpers",
"JSON::Any",
"File::Path",
"File::NFSLock",
"IO::Dir",
"Path::Class",
"Carp::Clan",
"MooseX::Types",
"MooseX::Types::Path::Class",
"Test::TempDir",
"Test::NoWarnings",
"Test::Tester",
"Test::Deep",
"String::RewritePrefix",
"MooseX::Storage",
"Math::Complex",
"Check::ISA",
"Hash::Util::FieldHash::Compat",
"Algorithm::C3",
"Class::C3",
"MRO::Compat",
"Scope::Guard",
"Devel::GlobalDestruction",
"Class::MOP",
"Try::Tiny",
"Moose",
"Task::Weaken",
"XSLoader",
"base",
"Variable::Magic",
"Data::OptList",
"Sub::Install",
"Params::Util",
"Sub::Exporter",
"B::Hooks::EndOfScope",
"Sub::Identify",
"Sub::Name",
"Package::Stash",
"namespace::clean",
"Tie::RefHash",
"Test::use::ok",
"Tie::ToObject",
"Data::Visitor",
"MooseX::Clone",
"Geometry::Primitive",
"ExtUtils::MakeMaker",
"Layout::Manager",
"Time::Local",
"List::MoreUtils",
"DateTime::Locale",
"File::Temp",
"Exporter",
"ExtUtils::ParseXS",
"Module::Build",
"Attribute::Handlers",
"ExtUtils::CBuilder",
"Params::Validate",
"Class::Singleton",
"DateTime::TimeZone",
"Test",
"Text::Wrap",
"Pod::Escapes",
"Pod::Simple",
"File::Spec",
"Pod::Man",
"Sub::Uplevel",
"Test::Exception",
"Test::Harness",
"Test::More",
"Scalar::Util",
"DateTime",
"Chart::Clicker");
