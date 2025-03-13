find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_program(SPHINX_APIDOC sphinx-apidoc REQUIRED)
find_program(FORD_EXE ford
  DOC "Path to the ford executable (required to generate the documentation)")
find_program(SPHINX_BUILD sphinx-build
  DOC "Path to sphinx-build executable (for generating Python documentation)")

# Ford
set(FORD_PROJECT_FILE "${CMAKE_BINARY_DIR}/docs/ford/openqp.md")
configure_file(${CMAKE_SOURCE_DIR}/docs/ford/openqp.md "${FORD_PROJECT_FILE}")

# Sphinx source (static docs in the source tree)
set(SPHINX_SOURCE "${CMAKE_SOURCE_DIR}/docs/sphinx")

# Define a build directory for Sphinx sources (including the generated API docs)
set(SPHINX_BUILD_SOURCE_DIR "${CMAKE_BINARY_DIR}/docs/sphinx_src")
# Define the output directory for the final HTML docs
set(SPHINX_BUILD_DIR "${CMAKE_SOURCE_DIR}/doc/sphinx")

set(PYTHON_SOURCE_DIR "${CMAKE_SOURCE_DIR}/pyoqp/oqp")

# Step 1: Copy static Sphinx docs from the source tree to the build directory.
add_custom_command(
  OUTPUT ${SPHINX_BUILD_SOURCE_DIR}/index.rst
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${SPHINX_SOURCE} ${SPHINX_BUILD_SOURCE_DIR}
  DEPENDS ${SPHINX_SOURCE}
  COMMENT "Copying Sphinx source files to build directory"
)

# Step 2: Generate API docs (modules.rst) in the build directory.
add_custom_command(
  OUTPUT ${SPHINX_BUILD_SOURCE_DIR}/modules.rst
  COMMAND ${SPHINX_APIDOC}
    -o ${SPHINX_BUILD_SOURCE_DIR}
    ${PYTHON_SOURCE_DIR}
  DEPENDS ${PYTHON_SOURCE_DIR} ${SPHINX_BUILD_SOURCE_DIR}/index.rst
  COMMENT "Running sphinx-apidoc to generate Python API docs"
  VERBATIM
)

add_custom_target(apidoc
  DEPENDS ${SPHINX_BUILD_SOURCE_DIR}/modules.rst
)

add_custom_target(ford_doc
  COMMENT "Generating Fortran API documentation with FORD"
  COMMAND "${FORD_EXE}" "${FORD_PROJECT_FILE}"
  VERBATIM
)

# Step 3: Build the final Sphinx HTML documentation using the build source directory.
add_custom_target(sphinx_doc
  COMMENT "Generating Python API documentation with Sphinx"
  COMMAND ${Python3_EXECUTABLE} -m sphinx -b html "${SPHINX_BUILD_SOURCE_DIR}" "${SPHINX_BUILD_DIR}"
  DEPENDS apidoc
  VERBATIM
)

configure_file(
  ${CMAKE_SOURCE_DIR}/docs/index.html  # Source index.html
  ${CMAKE_SOURCE_DIR}/doc/index.html   # Destination (build output)
  COPYONLY                             # Just copy without parsing
)

add_custom_target(doc
  COMMENT "Generating complete API documentation (Python + Fortran)"
  DEPENDS ford_doc sphinx_doc
)

