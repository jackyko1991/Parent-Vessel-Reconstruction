list(APPEND PROJECT_LIBRARIES 
	Core
	)

foreach(ONE ${PROJECT_LIBRARIES})
	add_subdirectory(${ONE})
	include_directories(${${ONE}_SOURCE_DIR})
	include_directories(${${ONE}_BINARY_DIR})
endforeach()
