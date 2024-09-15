# Contributing

### Step-by-Step Guide

#### 1. **Fork the Repository**
Forking creates a personal copy of the MicroMagnetic.jl repository under your GitHub account.

- Go to the **MicroMagnetic.jl** repository: [MicroMagnetic.jl](https://github.com/magneticsimulation/MicroMagnetic.jl).
- In the top-right corner, click the **Fork** button to create a copy of the repository in your own GitHub account.

#### 2. **Clone Your Fork**
After forking the repository, you need to clone your fork to your local machine to start working on it.

- Open your terminal (or Git Bash on Windows).
- Clone the repository using the following command (replace `yourusername` with your GitHub username):
  ```bash
  git clone https://github.com/yourusername/MicroMagnetic.jl
  cd MicroMagnetic.jl
  ```

#### 3. **(Optional) Set the Original Repository as Upstream**
To keep your fork in sync with the original repository, you can add the original `MicroMagnetic.jl` repository as a remote named `upstream`. This step is optional, but it helps you easily pull in updates from the main repository.

- Run the following commands to set it up:
  ```bash
  git remote add upstream https://github.com/magneticsimulation/MicroMagnetic.jl
  git fetch upstream
  ```

If you don't set up `upstream`, you can still contribute, but you’ll need to manually manage any updates from the original repository.

#### 4. **Create a New Branch**
Before making any changes, create a new branch to isolate your modifications from the main branch. It's good practice to create descriptive branch names, like `fix-simulation-bug` or `add-feature-x`.

- Create and switch to a new branch:
  ```bash
  git checkout -b my-new-feature
  ```

#### 5. **Activate Local Development in Julia**
To ensure that Julia uses your local version of the package rather than the registered one, activate development mode in Julia. Assuming you have cloned the repository to `/path/to/MicroMagnetic.jl`, use the following command:

```julia
  using Pkg
  Pkg.develop(path="/path/to/MicroMagnetic.jl")
```

This allows you to work on the local version of the package.

#### 6. **Make Your Changes**
Now that you are in your new branch, you can start making changes to the code. For example, you can add new functionality, fix bugs, or update documentation.

- Use your favorite text editor or IDE to modify the source code in the **MicroMagnetic.jl** repository.

#### 7. **Test Your Changes**
Before pushing your changes, ensure they work correctly by running tests. The tests are usually located in the `test` folder.

- To run tests:
  ```bash
  julia --project test/runtests.jl
  ```

Make sure all tests pass before proceeding. If any tests fail, debug and fix your code.

#### 8. **Commit Your Changes**
Once you're happy with the changes, you need to commit them. Make sure to write a clear and concise commit message describing what you’ve done.

- Stage your changes:
  ```bash
  git add .
  ```
- Commit with a message:
  ```bash
  git commit -m "Add feature X to improve simulation accuracy"
  ```

#### 9. **Push Your Changes to Your Fork**
Now that you’ve committed your changes locally, push them to your forked repository on GitHub.

- Push your branch:
  ```bash
  git push origin my-new-feature
  ```

#### 10. **Open a Pull Request**
Once your changes are pushed to your fork on GitHub, you’re ready to open a pull request (PR) to merge your changes back into the main repository.

- Go to your forked repository on GitHub.
- You’ll see a button prompting you to open a pull request. Click on it.
- In the pull request, provide a detailed description of the changes you made, why they are needed, and any relevant information.
- Submit the pull request.

#### 11. **(Optional) Keep Your Fork Up-to-Date**
As the original `MicroMagnetic.jl` repository evolves, you’ll want to keep your fork in sync with the latest changes. If you followed step 3 to set `upstream`, you can do this easily.

- Fetch the latest updates from the original repository:
  ```bash
  git fetch upstream
  ```
- Merge the changes into your local branch:
  ```bash
  git merge upstream/main
  ```

### Summary of Key Commands
Here’s a quick reference for the key Git commands used in this guide:

- Fork: (via GitHub UI)
- Clone: `git clone https://github.com/yourusername/MicroMagnetic.jl`
- (Optional) Set Upstream: `git remote add upstream https://github.com/magneticsimulation/MicroMagnetic.jl`
- New Branch: `git checkout -b my-new-feature`
- Add Changes: `git add .`
- Commit: `git commit -m "message"`
- Push: `git push origin my-new-feature`
- (Optional) Fetch Updates: `git fetch upstream`
- (Optional) Merge Updates: `git merge upstream/main`
- Activate Local Development: `Pkg.develop(path="/path/to/MicroMagnetic.jl")`

Feel free to explore, experiment, and contribute to making **MicroMagnetic.jl** even better!