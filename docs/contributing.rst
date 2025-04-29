============
Contributing
============

We welcome your input and contributions to pyROX! If you have suggestions for new features, or if you have found a (potential) bug, do not hesitate to let us know by `opening an issue on GitHub <https://github.com/samderegt/pyROX/issues>`_.

If you would like to implement a new feature or fix a bug, please follow these steps:
1. **Fork the repository**: Click the "Fork" button at the top right of the `GitHub page <https://github.com/samderegt/pyROX>`_ to create your own copy of the repository.
2. **Install the repository**: Install your forked repository to your machine:
   .. code-block:: bash

      git clone https://github.com/<your-user-name>/pyROX
      pip install -e .

   This will install the package in editable mode, allowing you to make changes to the source code and see them reflected in your Python environment immediately.
3. **Create a new branch**: Before making changes, it is best to create a new branch:
   .. code-block:: bash

      git checkout -b feature/your-feature-name
   
4. **Make your changes**: Implement the new feature and push to your forked repository.
5. **Create a pull request**: Go to the new feature's branch on GitHub (in your forked repository) and click "Contribute" > "Open pull request". Select ``samderegt/pyROX:dev`` as the base branch to merge into, and your feature branch as the head branch. Check and resolve any conflicts and add a description of your changes. When you are ready, click "Create pull request".
   .. note:: We will implement new features in the `dev` branch. Once a feature is stable, we will make a new release and merge the changes into the `main` branch.
   .. note:: Creating the pull request will automatically trigger some tests to ensure that the changes do not break any existing functionality. You can also run the tests locally by executing ``pytest`` in the root directory.
6. **Review and discuss**: We will now review your pull request and may ask you to make changes or provide additional information. Once the changes are satisfactory, we will merge the pull request into the ``dev`` branch. Thank you for your contribution to pyROX!