<!--
  This directory is the canonical source for the `lean/` copy inside the
  P1, P2 and P4 packets: `packets.py build` mirrors every file here into all
  three and the gate pins their SHA-256. Editing any file in `lean/` --- this
  README included, boilerplate though it looks --- turns `P1` red on three
  packets at once until `packets.py build` is run. Edit here, never in a
  packet's copy, and rebuild in the same change.
-->

# goldbach-lean

## GitHub configuration

To set up your new GitHub repository, follow these steps:

* Under your repository name, click **Settings**.
* In the **Actions** section of the sidebar, click "General".
* Check the box **Allow GitHub Actions to create and approve pull requests**.
* Click the **Pages** section of the settings sidebar.
* In the **Source** dropdown menu, select "GitHub Actions".

After following the steps above, you can remove this section from the README file.
