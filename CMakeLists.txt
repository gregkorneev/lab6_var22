cmake_minimum_required(VERSION 3.16)

project(gauss_jordan LANGUAGES CXX)

# C++17 без расширений компилятора
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Тип сборки по умолчанию: Release (если не задано извне)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Исполняемый файл
add_executable(app src/main.cpp)

# Кладём бинарники в build/bin
set_target_properties(app PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Предупреждения
if(MSVC)
  target_compile_options(app PRIVATE /W4 /permissive-)
else()
  target_compile_options(app PRIVATE -Wall -Wextra -Wpedantic)
endif()

# LTO/IPO если поддерживается
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_msg)
if(ipo_supported)
  set_property(TARGET app PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Установка (опционально: make install / cmake --install)
install(TARGETS app RUNTIME DESTINATION bin)
