add_executable(tests main.cc)

target_link_libraries(tests PUBLIC Catch2::Catch2WithMain fmt::fmt httplib::httplib)
target_compile_features(tests PUBLIC cxx_std_20)

catch_discover_tests(tests)
