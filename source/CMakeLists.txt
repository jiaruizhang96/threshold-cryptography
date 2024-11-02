add_executable(main main.cc)

target_link_libraries(main PUBLIC fmt::fmt httplib::httplib)
target_compile_features(main PUBLIC cxx_std_20)
