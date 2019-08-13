#include <iostream>
#include <protocols/features/FeaturesReporterFactory.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <core/init/init.hh>
#include <vector>
#include <string>
#include <devel/init.hh>
#include <core/types.hh>
#include <utility/excn/Exceptions.hh>

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols;
		using namespace protocols::features;
		devel::init(argc,argv);
		protocols::features::FeaturesReporterFactory *factory = FeaturesReporterFactory::get_instance();
		utility::vector1<std::string> all_features = factory->get_all_features_names();
		for ( core::Size i=1; i<=all_features.size(); ++i ) {
			FeaturesReporterOP feature_reporter = factory->get_features_reporter(all_features[i]);
			std::string schema = feature_reporter->schema();
			std::cout<<schema<<std::endl;
		}

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
