#!/bin/bash

echo "GVClass Pipeline Threading Benchmarks (with pixi)"
echo "=================================================="
echo ""

# Run 1: IQ-TREE with 32 threads
echo "Test 1: IQ-TREE with 32 threads"
echo "--------------------------------"
rm -rf example_results/
START=$(date +%s)
pixi run ./gvclass example --tree-method iqtree -t 32
END=$(date +%s)
RUNTIME1=$((END-START))
echo "Runtime: ${RUNTIME1} seconds"
echo ""

# Run 2: IQ-TREE with 8 threads  
echo "Test 2: IQ-TREE with 8 threads"
echo "-------------------------------"
rm -rf example_results/
START=$(date +%s)
pixi run ./gvclass example --tree-method iqtree -t 8
END=$(date +%s)
RUNTIME2=$((END-START))
echo "Runtime: ${RUNTIME2} seconds"
echo ""

# Run 3: FastTree with 32 threads
echo "Test 3: FastTree with 32 threads"
echo "---------------------------------"
rm -rf example_results/
START=$(date +%s)
pixi run ./gvclass example --tree-method fasttree -t 32
END=$(date +%s)
RUNTIME3=$((END-START))
echo "Runtime: ${RUNTIME3} seconds"
echo ""

# Run 4: FastTree with 8 threads
echo "Test 4: FastTree with 8 threads"
echo "--------------------------------"
rm -rf example_results/
START=$(date +%s)
pixi run ./gvclass example --tree-method fasttree -t 8
END=$(date +%s)
RUNTIME4=$((END-START))
echo "Runtime: ${RUNTIME4} seconds"
echo ""

echo "=================================================="
echo "Summary of Results:"
echo "=================================================="
echo "IQ-TREE  32 threads: ${RUNTIME1} seconds"
echo "IQ-TREE   8 threads: ${RUNTIME2} seconds"
echo "FastTree 32 threads: ${RUNTIME3} seconds"
echo "FastTree  8 threads: ${RUNTIME4} seconds"
echo ""

# Calculate speedups
if [ $RUNTIME2 -gt 0 ]; then
    SPEEDUP_IQ=$(echo "scale=2; $RUNTIME2 / $RUNTIME1" | bc -l 2>/dev/null || echo "N/A")
    echo "IQ-TREE speedup (32 vs 8 threads): ${SPEEDUP_IQ}x"
fi

if [ $RUNTIME4 -gt 0 ]; then
    SPEEDUP_FT=$(echo "scale=2; $RUNTIME4 / $RUNTIME3" | bc -l 2>/dev/null || echo "N/A")
    echo "FastTree speedup (32 vs 8 threads): ${SPEEDUP_FT}x"
fi