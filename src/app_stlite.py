from varlociraptor_inspect.views.main import main_view

await main_view()  # noqa: F704, PLE1142 — stlite compiles with top-level await support
