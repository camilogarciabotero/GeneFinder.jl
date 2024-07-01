@testset "Aqua tests" begin

    Aqua.test_ambiguities(GeneFinder)
    Aqua.test_persistent_tasks(GeneFinder)
    Aqua.test_piracies(GeneFinder)
    Aqua.test_stale_deps(GeneFinder)
    # Aqua.test_unbound_args(GeneFinder)
    Aqua.test_undefined_exports(GeneFinder)

end