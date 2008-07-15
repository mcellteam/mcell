from test_parser import oldvizsuite, vizsuite, errorsuite, quicksuite, kitchensinksuite, rtcheckpointsuite

tests = {
    "oldvizsuite"       : "VIZ output tests for ASCII/RK/DX modes",
    "vizsuite"          : "VIZ output tests for DREAMM V3 modes",
    "errorsuite"        : "Test error handling for invalid MDL files",
    "quicksuite"        : "A few quick running tests which cover most valid MDL options",
    "kitchensinksuite"  : "Kitchen Sink Test: (very nearly) every parser option",
    "rtcheckpointsuite" : "Basic test of timed checkpoint functionality"
}
collections = {
    "allvizsuite"  : ("VIZ output tests for all modes (old+new)", ["oldvizsuite", "vizsuite"]),
    "fasttests"    : ("All quick running tests (valid+invalid MDL)", ["errorsuite", "quicksuite"]),
}
