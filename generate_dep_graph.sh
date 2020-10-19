#!/bin/bash

rm -rf build_debug_for_dep_graph 
mkdir -p build_debug_for_dep_graph

pushd build_debug_for_dep_graph && cmake -DCMAKE_BUILD_TYPE=Debug --graphviz=fig-dep-all.dot .. && popd
pushd build_debug_for_dep_graph && dot -Tpng fig-dep-all.dot.monostar_app -o ../fig-dep-all.png && popd

echo 'set(GRAPHVIZ_EXTERNAL_LIBS false)' > CMakeGraphVizOptions.cmake

pushd build_debug_for_dep_graph && cmake -DCMAKE_BUILD_TYPE=Debug --graphviz=fig-dep-no-enxternal.dot .. && popd
pushd build_debug_for_dep_graph && dot -Tpng fig-dep-no-enxternal.dot.monostar_app -o ../fig-dep-no-enxternal.png && popd

rm -rf build_debug_for_dep_graph 

