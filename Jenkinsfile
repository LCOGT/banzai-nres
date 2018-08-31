#!/usr/bin/env groovy

@Library('lco-shared-libs@0.0.8') _

pipeline {
	agent any
	environment {
		dockerImage = null
		PROJ_NAME = projName()
		GIT_DESCRIPTION = gitDescription()
		DOCKER_IMG = dockerImageName("${LCO_DOCK_REG}", "${PROJ_NAME}", "${GIT_DESCRIPTION}")
	}
	options {
		timeout(time: 3, unit: 'HOURS')
	}
	stages {
		stage('Build image') {
			steps {
				script {
					dockerImage = docker.build("${DOCKER_IMG}", "--build-arg BANZAI_VERSION=0.9.7-60-gd81e26b .")
				}
			}
		}
		stage('Push image') {
			steps {
				script {
					dockerImage.push("${GIT_DESCRIPTION}")
				}
			}
		}
		stage('Test') {
			steps {
				script {
					sh 'docker run --rm --user=root --entrypoint=pytest ${DOCKER_IMG} -m "not e2e" /lco/banzai-nres/'
				}
			}
		}
		stage('DeployTestStack') {
			when {
				anyOf {
					branch 'PR-*'
					expression { return params.forceEndToEnd }
				}
			}
			environment {
				RANCHERDEV_CREDS = credentials('rancher-cli-dev')
			}
			steps {
				script {
					withCredentials([usernamePassword(
							credentialsId: 'rabbit-mq',
							usernameVariable: 'RABBITMQ_USER',
							passwordVariable: 'RABBITMQ_PASSWORD')]) {
						sh('rancher -c ${RANCHERDEV_CREDS} up --stack BANZAINRESPipelineTest --force-upgrade --confirm-upgrade -d')
					}
				}
			}
		}
		stage('TestE2E') {
			when {
				anyOf {
					branch 'PR-*'
					expression { return params.forceEndToEnd }
				}
			}
			environment {
				RANCHERDEV_CREDS = credentials('rancher-cli-dev')
				SSH_CREDS = credentials('jenkins-rancher-ssh-userpass')
				CONTAINER_ID = getContainerId('BANZAINRESPipelineTest-BANZAINRESPipelineTest-1')
				CONTAINER_HOST = getContainerHostName('BANZAINRESPipelineTest-BANZAINRESPipelineTest-1')
				ARCHIVE_UID = credentials('archive-userid')
			}
			steps {
				script {
					sshagent(credentials: ['jenkins-rancher-ssh']) {
						executeOnRancher('pytest -m e2e /lco/banzai-nres/', CONTAINER_HOST, CONTAINER_ID, ARCHIVE_UID)
					}
				}
			}
			post {
				success {
					script {
						sh('rancher -c ${RANCHERDEV_CREDS} rm --stop --type stack BANZAINRESPipelineTest ')
					}
				}
			}
		}
	}
}
