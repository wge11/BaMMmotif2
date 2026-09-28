# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------

# command-line parsing (third-party GetOpt_pp, GPLv3)
add_library(getopt_pp STATIC getopt_pp/getopt_pp.cpp)
target_include_directories(getopt_pp PUBLIC getopt_pp)

# alphabet, sequences, background models and motif models
add_library(bamm_init STATIC
    init/Alphabet.cpp
    init/BackgroundModel.cpp
    init/BackgroundModelSet.cpp
    init/BackgroundModelSetScore.cpp
    init/Motif.cpp
    init/MotifSet.cpp
    init/Sequence.cpp
    init/SequenceSet.cpp
)
target_include_directories(bamm_init PUBLIC init refinement)

# model optimisation (EM, Gibbs sampling), scoring, negative-set sampling
# and cross-validation shared by all command-line tools
add_library(bamm_core STATIC
    refinement/EM.cpp
    refinement/GibbsSampling.cpp
    seq_generator/SeqGenerator.cpp
    seq_scoring/ScoreSeqSet.cpp
    evaluation/FDR.cpp
)
target_link_libraries(bamm_core PUBLIC bamm_init)

# ---------------------------------------------------------------------------
# Executables
# ---------------------------------------------------------------------------

add_executable(BaMMmotif refinement/mainBaMM.cpp refinement/Global.cpp)
target_link_libraries(BaMMmotif PRIVATE bamm_core getopt_pp)

add_executable(BaMMScan seq_scoring/mainScan.cpp seq_scoring/GScan.cpp)
target_link_libraries(BaMMScan PRIVATE bamm_core)

add_executable(FDR evaluation/mainFDR.cpp evaluation/GFdr.cpp)
target_link_libraries(FDR PRIVATE bamm_core)

add_executable(BaMMSimu seq_generator/mainSimulator.cpp seq_generator/GSimu.cpp)
target_link_libraries(BaMMSimu PRIVATE bamm_core)

add_executable(extractProbs misc/extractProbs.cpp)
target_link_libraries(extractProbs PRIVATE bamm_init)

install(TARGETS BaMMmotif BaMMScan FDR BaMMSimu extractProbs
        RUNTIME DESTINATION bin)
