# CurveSmash website

The static website source lives here, alongside the challenge definitions in
`Layr-Labs/curvesma.sh`. It shows all four tracks, their current cost rules and
the shared verifier's submission process. No scores are published until reviewed
under the current rules. Historical evidence remains linked separately.

Run `python3 scripts/build_website.py` to export the protected profiles and cost
model into `website/dist/`, then preview:

```sh
python3 -m http.server 4173 --bind 127.0.0.1 --directory website/dist
```

The Pages workflow builds and publishes this directory on changes to `master`.
GitHub Pages must first be enabled with **GitHub Actions** as the source. The site must be **public**, including when the repository is internal. A
repository administrator must move the custom domain `curvesma.sh` from
`zac-williamson/curvesmash` to this repository and enable HTTPS. Keep the existing
site running until the replacement is ready. Both repositories belong to GitHub
Pages, so the existing apex DNS addresses can remain unchanged.

`dist/judges.json` lists appointed human judges, each with `name` and
`affiliation`. An empty or unavailable roster hides the banner. This display
does not grant reviewer access in Yukon. Registration and opening the four
tracks in Yukon remain separate from website hosting.
