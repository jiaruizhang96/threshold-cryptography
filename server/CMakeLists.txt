add_executable(server main.cc)

target_link_libraries(server PUBLIC fmt::fmt httplib::httplib nlohmann_json::nlohmann_json)
target_compile_features(server PUBLIC cxx_std_20)
