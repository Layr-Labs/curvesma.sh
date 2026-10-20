"""Keep the public presentation aligned with protected challenge definitions."""
from html.parser import HTMLParser
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / 'website' / 'dist'


class Page(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.elements = []
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        self.elements.append((tag, dict(attrs)))


class Website(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('build_website',
                                                     ROOT / 'scripts/build_website.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.build()

    def test_all_tracks_have_cards_selectors_and_submission_buttons(self):
        tracks = json.loads((ROOT / 'benchmark.json').read_text())['tracks']
        names = {track['name'] for track in tracks}
        page = Page((SITE / 'index.html').read_text())
        self.assertEqual({a['data-target'] for _, a in page.elements
                          if 'data-target' in a}, names)
        self.assertEqual({a['data-prepare'] for _, a in page.elements
                          if 'data-prepare' in a}, names)
        # Each track is selectable in both results and participation.
        for name in names:
            self.assertEqual(sum(tag == 'option' and a.get('value') == name
                                 for tag, a in page.elements), 2)
        roster = json.loads((SITE / 'judges.json').read_text())
        self.assertIsInstance(roster, list)
        for judge in roster:
            self.assertTrue(judge['name'].strip())
            self.assertIsInstance(judge['affiliation'], str)

    def test_public_definitions_are_exported_from_the_protected_sources(self):
        manifest = json.loads((ROOT / 'benchmark.json').read_text())
        for track in manifest['tracks']:
            name = track['name']
            self.assertEqual((SITE / 'targets' / f'{name}.json').read_bytes(),
                             (ROOT / 'tracks' / name / 'target-profile.json').read_bytes())
        model = json.loads((SITE / 'cost-model.json').read_text())
        self.assertEqual(model,
                         json.loads((ROOT / 'tracks/classical-group-work-v2.json').read_text()))
        html = (SITE / 'index.html').read_text()
        for price in set(model['reference_operation_costs'].values()):
            self.assertIn(f'<td>{price:,}</td>', html)
        self.assertIn('W(2²⁰)', html)
        self.assertIn('batched incomplete affine 6M', html)

    def test_local_links_resolve_and_obsolete_participation_is_removed(self):
        html = (SITE / 'index.html').read_text()
        js = (SITE / 'app.js').read_text()
        page = Page(html)
        ids = {a['id'] for _, a in page.elements if 'id' in a}
        for _, attrs in page.elements:
            for key in ('href', 'src'):
                url = attrs.get(key, '')
                if url.startswith('./'):
                    self.assertTrue((SITE / url[2:]).is_file(), url)
                if url.startswith('#'):
                    self.assertIn(url[1:], ids)
            if 'data-open' in attrs:
                self.assertIn(attrs['data-open'], ids)
        for obsolete in ('candidates/', 'local_tracks.py', '.yukon/setup.sh',
                         'Four AI reviews', 'gpt-5.6-sol', '129.000802', '128.494762'):
            self.assertNotIn(obsolete, html + js)
        self.assertIn('https://github.com/Layr-Labs/curvesma.sh.git', js)
        self.assertIn('selected-submission.txt', js)
        self.assertIn('score:{numerator:', js)


if __name__ == '__main__':
    unittest.main()
